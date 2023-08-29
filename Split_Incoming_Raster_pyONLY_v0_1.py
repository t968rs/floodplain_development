import sys

import arcpy
import time
import os

import concurrent.futures


# from PBAR_formatting import PBARformatting
# from tqdm import tqdm


class RasterSplitter:

    def __init__(self, wse_grid_folder, output_folder, tile_wse, split_by, huc12_path, watch):

        self.input_raster_grid_folder = wse_grid_folder
        self.output_folder = output_folder
        self.tile_input_rasters = tile_wse
        self.timer = watch
        self.huc12_polygon_path = huc12_path
        self.splitby = split_by

        self.cpu_count = int(round(64, 0))

        self.shrunken_raster_folder = ''
        self.scratch_gdb = ''
        self.system_memory = ''  # GB integer

        self.raster_grid_dictionary = {}
        self.shrunken_raster_dictionary = {}
        self.huc12_raster_dictionary = {}

        self.contour_fclist = []

        self.variable_type_conversions()

        # self.pbar_format = PBARformatting.get_barformat()

    def variable_type_conversions(self):

        # String -> Boolean
        if type(self.tile_input_rasters) == str:
            if self.tile_input_rasters == '' or self.tile_input_rasters == '#':
                self.tile_input_rasters = False
            else:
                # arcpy.AddMessage(f'Tile? {self.tile_input_rasters}')
                self.tile_input_rasters = True

        # Paths from inputs
        for var_name in ['huc12_polygon_path']:
            var_value = getattr(self, var_name)
            arcpy.AddMessage(f'# Var {var_name}: {var_value}')
            if var_value is not None:
                arcpy.AddMessage(f'  is not <Null>')
                if not os.path.split(var_value)[0]:
                    desc = arcpy.Describe(var_value)
                    base_path, name = desc.path, desc.name
                    new_path = os.path.join(base_path, name)
                    setattr(self, var_name, new_path)
                    arcpy.AddMessage(f'  Updated path for: {var_value}: {new_path}')

    def create_output_folder(self):

        if not os.path.exists(self.output_folder):
            os.makedirs(self.output_folder)

        # Create output folder
        shrunken_rasters_folder = os.path.join(self.output_folder, '01_tiled_rasters')
        if not os.path.exists(shrunken_rasters_folder):
            os.makedirs(shrunken_rasters_folder)
        self.shrunken_raster_folder = shrunken_rasters_folder

        # Create scratch GDB
        scratch_gdb = os.path.join(self.output_folder, 'scratch.gdb')
        if not arcpy.Exists(scratch_gdb):
            arcpy.management.CreateFileGDB(self.output_folder, 'scratch.gdb')

        self.scratch_gdb = scratch_gdb

    def create_raster_grid_list(self, folder):

        raster_grids = {}
        tiff_list = os.walk(folder)
        for root, folders, files in tiff_list:
            for file in files:
                if '.tif' in file.lower():
                    split = file.split('.tif')
                    # arcpy.AddMessage(f'File {file}: {split[1]}')
                    if split[1] == '':
                        tifname = file.split(".tif")[0]
                        raster_grids[tifname] = os.path.join(root, os.path.join(root, file))
        self.raster_grid_dictionary = raster_grids

    def populate_raster_props(self):

        new_tiff_dict = {}
        for tifname, tifpath in self.raster_grid_dictionary.items():
            arcpy.AddMessage(f'Raster: {tifpath}')
            file_size = round(os.stat(tifpath).st_size * 0.000001, 0)
            arcpy.AddMessage(f' {tifname} is {file_size} MB')
            desc_tiff = arcpy.Describe(tifpath)
            tiff_sr_wkid = desc_tiff.spatialReference.factoryCode
            cellsize = desc_tiff.meanCellHeight
            tiff_width = desc_tiff.width
            tiff_height = desc_tiff.height
            tiff_extent = desc_tiff.extent

            if tiff_width / tiff_height > 2 or tiff_height / tiff_width > 2:
                wh_multiplier = tiff_width / tiff_height
            else:
                wh_multiplier = None

            new_tiff_dict[tifname] = {'file mb': file_size, 'sr_wkid': tiff_sr_wkid,
                                      'cell size': cellsize, 'width': tiff_width, 'height': tiff_height,
                                      'path': tifpath, 'wh multiplier': wh_multiplier, 'extent': tiff_extent}
            arcpy.AddMessage(f' {new_tiff_dict}')

        self.raster_grid_dictionary = new_tiff_dict

    def convert_rasters(self):

        # Check raster properties and export as necessary (file size)
        ## Get pbar too!
        arcpy.env.overwriteOutput = True
        # pbar = tqdm(total=1 * len(self.wse_grid_dictionary), bar_format=self.pbar_format,
        # desc=f'Copying/Splitting WSE grids...')
        converted_count = 0
        converted_into = {}
        new_wse_dict = {}
        for tiffname, tiff_dictionary in self.raster_grid_dictionary.items():
            current_size_mb = tiff_dictionary['file mb']
            input_path = tiff_dictionary['path']
            wh_multiplier = tiff_dictionary['wh multiplier']
            number_tiles = int(current_size_mb / 20) + 1
            if (number_tiles % 2) != 0:
                number_tiles += 1

            if wh_multiplier is not None:
                tiles_x = int(round(number_tiles * wh_multiplier / 2))
                tiles_y = int(round(number_tiles / wh_multiplier / 2))
            else:
                tiles_x = int(round(number_tiles / 2))
                tiles_y = int(round(number_tiles / 2))

            if current_size_mb > 100 and tiles_y != 1 and tiles_x != 1:

                # pbar.set_postfix(splitting=f'New Tiles: {tiles_x}x{tiles_y} gridded tiles, {number_tiles} tiles')
                converted_into[tiffname] = number_tiles

                arcpy.management.SplitRaster(input_path, self.shrunken_raster_folder, out_base_name=tiffname,
                                             split_method='NUMBER_OF_TILES', format='TIFF', resampling_type='BILINEAR',
                                             num_rasters=f"{tiles_x} {tiles_y}", overlap=2, nodata_value=-9999)
            else:
                output_raster = os.path.join(self.shrunken_raster_folder, tiffname + "_.tif")
                arcpy.management.CopyRaster(input_path, output_raster,
                                            format='TIFF', pixel_type='32_BIT_SIGNED', nodata_value=-9999)

            # Add items to dictionary by searching for rasters
            arcpy.env.workspace = self.shrunken_raster_folder
            dictionary_value = self.raster_grid_dictionary[tiffname]
            tiff_list = arcpy.ListRasters(wild_card=f'*{tiffname}*')
            for raster in tiff_list:
                filename = os.path.split(raster)[1]
                name = filename.split('.')[0]
                new_wse_dict[name] = dictionary_value

            # pbar.update(1)
            converted_count += 1
        # del pbar
        self.raster_grid_dictionary.update(new_wse_dict)

        return converted_count, converted_into

    def create_huc12_dicts(self):

        # Check raster properties and export as necessary (file size)
        ## Get pbar too!
        arcpy.env.overwriteOutput = True
        # pbar = tqdm(total=1 * len(self.wse_grid_dictionary), bar_format=self.pbar_format,
        # desc=f'Copying/Splitting WSE grids...')
        huc12_raster_dict = {}

        # Find HUC12_ID and HUC12 Name fields
        huc12_fnames = [f.name for f in arcpy.ListFields(self.huc12_polygon_path)]
        huc_fname = ''
        name_fname = ''
        for fname in huc12_fnames:
            if 'huc12' in fname.lower():
                huc_fname = fname
            elif 'name' in fname.lower():
                name_fname = fname

        # Iterate incoming input rasters from dictionary

        for tiffname in list(self.raster_grid_dictionary.keys()):
            current_size_mb = self.raster_grid_dictionary[tiffname]['file mb']
            raster_path = self.raster_grid_dictionary[tiffname]['path']
            input_extent = self.raster_grid_dictionary[tiffname]['extent']

            scratch_huc12_fc = os.path.join(self.scratch_gdb, f'huc12_{tiffname}')

            if current_size_mb > 400:
                # Create subset of HUC12 that overlap the input dataset
                arcpy.management.CreateFeatureclass(self.scratch_gdb, f'huc12_{tiffname}',
                                                    template=self.huc12_polygon_path,
                                                    spatial_reference=self.huc12_polygon_path)
                huc12_inserts = arcpy.da.InsertCursor(scratch_huc12_fc,
                                                      ['Shape@', huc_fname, name_fname])

                huc12count = 0
                huc12_list = []
                with arcpy.da.SearchCursor(self.huc12_polygon_path, ['Shape@', huc_fname, name_fname]) as cursor:
                    for row in cursor:
                        if not row[0].disjoint(input_extent):
                            huc12_inserts.insertRow(row)
                            huc12_list.append(row[1])
                            huc12count += 1

                huc12_raster_dict[tiffname] = {'HUC12 List': huc12_list, 'HUC12 Temp Path': scratch_huc12_fc,
                                               'HUC 12 Field Name': huc_fname, 'Raster Path': raster_path}
                del huc12_inserts
                del self.raster_grid_dictionary[tiffname]['extent']

            else:  # If smaller than 400MB, just export to the new folder
                outrasterr_path = os.path.join(self.shrunken_raster_folder, tiffname)
                arcpy.management.CopyRaster(raster_path, outrasterr_path,
                                            format='TIFF', pixel_type='32_BIT_SIGNED', nodata_value=-9999)

        self.huc12_raster_dictionary = huc12_raster_dict

    def clip_rasters_manager(self):
        output_dictionary = {}
        processing_times = {}
        # Send concurrent stuff to the processor using many processes
        for tiff, tiff_dictionary in self.huc12_raster_dictionary.items():
            huc12_list = self.huc12_raster_dictionary[tiff]['HUC12 List']
            huc12_temppath = self.huc12_raster_dictionary[tiff]['HUC12 Temp Path']
            huc_fname = self.huc12_raster_dictionary[tiff]['HUC 12 Field Name']
            raster_path = self.huc12_raster_dictionary[tiff]['Raster Path']

            n_workers = min([int(round(self.cpu_count * 0.3, 0)), len(huc12_list)])
            chunksize = 1
            arcpy.AddMessage(
                f"  Using {n_workers} logical processors to clip rasters, {chunksize} HUC12 at a time")
            arcpy.AddMessage(f'    Passing {len(huc12_list)} HUC IDs to parallel processor:')
            arcpy.AddMessage(f'    {huc12_list}')
            arcpy.AddMessage(f'    ...for HUC12 watersheds in {huc12_temppath}')
            outras_pathlist = []
            complete_huc12_list = []
            with concurrent.futures.ProcessPoolExecutor(n_workers) as pool:
                futurelist = []
                for huc12 in huc12_list:

                    futurelist.append(pool.submit(self.raster_clip_worker, folder=self.shrunken_raster_folder,
                                                  pg_id=huc12, polygons=huc12_temppath,
                                                  fname=huc_fname, raster_path=raster_path))
                    processing_times = {f'Start Clip {huc12}': time.time()}
                for future in concurrent.futures.as_completed(futurelist, timeout=300):
                    outraster_path, huc_id_complete = future.result()
                    arcpy.AddMessage(f'  Finished {os.path.split(raster_path)[1]} raster export')
                    outras_pathlist.append(outraster_path)
                    complete_huc12_list.append(huc_id_complete)
                    processing_times[f'Finish Clip {huc12}'] = time.time()
            output_raster_info = (complete_huc12_list, outras_pathlist)

            output_dictionary = {tiff: {'Shrunken Raster Tuple': output_raster_info}}

        self.shrunken_raster_dictionary.update(output_dictionary)
        return processing_times

    @staticmethod
    def raster_clip_worker(folder, pg_id, polygons, fname, raster_path):

        # Make temp copy of polygon layer
        temp_pg = os.path.join('memory', f'pg_{pg_id}')
        if arcpy.Exists(temp_pg):
            arcpy.management.Delete(temp_pg)
        sql = f"{fname} = '{pg_id}'"
        mem_copy = arcpy.conversion.ExportFeatures(polygons, temp_pg, where_clause=sql)
        geo = arcpy.management.CopyFeatures(mem_copy, arcpy.Geometry())[0]

        # Clip raster by HUC12 polygons
        raster_name = f'h12_{str(pg_id)}.tif'
        rectangle = geo.hullRectangle
        out_path = os.path.join(folder, raster_name)
        if not os.path.exists(out_path):
            # clipping_geo = arcpy.management.CopyFeatures(row[0], arcpy.Geometry)
            # hucfname_formatted = arcpy.AddFieldDelimiters(mem_copy, huc_fname)
            # arcpy.AddMessage(f' SQL: <{sql}>')
            if arcpy.Exists('huc12_temp'):
                arcpy.management.Delete('huc12_temp')
            temp_huc12 = arcpy.management.MakeFeatureLayer(mem_copy, 'huc12_temp',
                                                           where_clause=sql)
            arcpy.AddMessage(f' Clipping {os.path.split(raster_path)[1]} by {pg_id} and outputing:\n   {out_path}')
            arcpy.management.Clip(
                in_raster=raster_path,
                rectangle=rectangle,
                out_raster=out_path,
                in_template_dataset=temp_huc12,
                nodata_value="-9999",
                clipping_geometry="ClippingGeometry",
                maintain_clipping_extent="NO_MAINTAIN_EXTENT"
            )
        else:
            arcpy.AddMessage(f'-- {out_path} alreayd exists')

        return str(out_path), str(pg_id)

    def confirm_output_rasters(self):

        # Get HUC12 names from output raster file names and append to list
        arcpy.env.workspace = self.shrunken_raster_folder
        output_foundraster_paths = [os.path.join(self.shrunken_raster_folder, r) for r in arcpy.ListRasters()]
        found_huc12_ids = []
        for outpath in output_foundraster_paths:
            name1 = outpath.split(".tif")[0]
            huc12name = name1.split("h12_")[1]
            found_huc12_ids.append(huc12name)

        # Check for expected HUC12s from input analysis
        for in_raster, related_info in self.huc12_raster_dictionary.items():
            huc12_expected_list = related_info['HUC12 List']
            huc12_confirmed_list = []
            for huc12_expected in huc12_expected_list:
                if huc12_expected not in found_huc12_ids:
                    arcpy.AddMessage(f' {huc12_expected} not found in output folder')
                else:
                    huc12_confirmed_list.append(huc12_expected)
            related_info['HUC12 Confirmed List'] = huc12_confirmed_list

        # Check to see if the number of clipped outputs, above, matches the number of expected clips
        for in_raster, related_info in self.huc12_raster_dictionary.items():
            confirmed_number = len(related_info['HUC12 Confirmed List'])
            expected_outputs = len(related_info['HUC12 List'])
            if confirmed_number >= expected_outputs:
                self.huc12_raster_dictionary[in_raster]['Confirmed in Output'] = True
            else:
                arcpy.AddMessage(f' Raster, {in_raster} confirmed {confirmed_number} outputs, but needed {expected_outputs}')

        # Print any rasters that did not create the expected number of clipped rasters.
        for raster, info in self.huc12_raster_dictionary.items():
            if not info['Confirmed in Output']:
                arcpy.AddMessage(f'Raster {raster} did not complete clipping by all overlappiung HUC12s')
                arcpy.AddMessage(f'  Info: {info}')

    def calc_total_raster_size(self):

        total_size = 0
        for name, properties in self.raster_grid_dictionary.items():
            size_mb = self.raster_grid_dictionary[name]['file mb']
            total_size = total_size + size_mb

        return total_size

    def run_raster_splitter(self):

        time_records = {'Start': time.time()}

        self.create_output_folder()
        self.create_raster_grid_list(folder=self.input_raster_grid_folder)
        arcpy.AddMessage(f'There are {len(self.raster_grid_dictionary)} WSE grids in the input folder')
        self.populate_raster_props()
        time_records['List/Properties Populated'] = time.time()
        total_raster_size = self.calc_total_raster_size()
        time_records[f'Total Raster Size {total_raster_size} MB'] = time.time()
        self.timer.time_reporter(times=time_records, new_iter=True)

        if self.splitby == 'Auto Tiles':
            count_converted, converted_dict = self.convert_rasters()
            time_records['Split All Rasters'] = time.time()
            arcpy.AddMessage(f'\nConverted {count_converted} rasters')
            for name, number in converted_dict.items():
                arcpy.AddMessage(f'{name} became {number} different rasters')
        else:
            self.create_huc12_dicts()
            raster_processing_times = self.clip_rasters_manager()
            time_records.update(raster_processing_times)
            time_records['Split All Rasters'] = time.time()
        self.timer.time_reporter(times=time_records, new_iter=False)

        self.confirm_output_rasters()
        arcpy.AddMessage(f'\nConverted {len(self.shrunken_raster_dictionary)} rasters')
        if self.splitby != 'Auto Tiles':
            for rastername, dictionary in self.shrunken_raster_dictionary.items():
                huc12s_overlap_raster = len(dictionary['Shrunken Raster Tuple'][0])
                huc12_output_confimed = self.huc12_raster_dictionary[rastername]['Confirmed in Output']
                arcpy.AddMessage(f'{rastername} overlaps {huc12s_overlap_raster} HUC12 polygons')
                if huc12_output_confimed:
                    arcpy.AddMessage(f'{rastername} successfully split')
                else:
                    arcpy.AddMessage(f'    ERROR: {rastername} did not complete all clips.')

        self.timer.time_reporter(times=time_records, new_iter=False)


class TimeReports:

    def __init__(self, workspace, options):

        self.times_dictionary = None
        self.start_time = None

        self.workspace = workspace
        self.options = options

    def get_start_time(self):

        start_time = 0
        for heading, a_time in self.times_dictionary.items():
            if 'start' in heading.lower():
                start_time = a_time
                break
        self.start_time = start_time

    def append_dictionary(self, times):

        if not self.times_dictionary:
            self.times_dictionary = times
        else:
            self.times_dictionary = self.times_dictionary.update(times)

    def time_reporter(self, times, new_iter):

        # Check for existing values
        self.append_dictionary(times)
        if not self.start_time:
            self.get_start_time()

        time_results = os.path.join(self.workspace, "time_results.txt")
        time_printout = open(time_results, "a")
        if new_iter:
            time_printout.writelines(f"\n\n    ---- THIS IS A NEW ITERATION ----\n\n")
            for option_def, option in self.options.items():
                time_printout.writelines(f'\n{option_def}: {option}')
        time_printout.close()

        timenames = list(times.keys())
        elapsed_times = {}
        self.times_dictionary = self.times_dictionary.update(times)

        for timename in timenames:
            total_seconds = times[timename]
            elapsed = total_seconds - self.start_time
            hours = elapsed // 3600
            elapsed = elapsed - 3600 * hours
            minutes = elapsed // 60
            elapsed = elapsed - 60 * minutes
            seconds = round(elapsed, 2)
            elapsed_times[timename] = (hours, minutes, seconds)
            print(f"{timename} processing finished after: {hours} hours, {minutes} minutes, {seconds} seconds")
            time_printout = open(time_results, "a")
            time_printout.writelines(f"\n{timename}: \n{hours} hours, {minutes} minutes, {seconds} seconds\n")
            time_printout.close()


if __name__ == '__main__':
    arcpy.AddMessage(' This script will fill and expand a raster, per user spec')

    input_raster_grid_folder = r'A:\2D_workflows\01_data\cochise'  # Paste a
    # folder location between quotes
    out_folder = r'A:\2D_workflows\03_test_outputs\cochise'  # Paste a folder location between quotes sys.argv[2]  #
    tile_input_wse = True
    split_raster_by = 'HUC12'  # sys.argv[3]  #
    split_huc12_path = r'A:\carto\HUC\AZ\WBDHU12.shp'  # sys.argv[4]  #

    run_options_m = {'Input Grid Folder': input_raster_grid_folder,
                     'Results Folder': out_folder,
                     'Generate Tiles?': tile_input_wse,
                     'Split By': split_raster_by,
                     'Split Path': split_huc12_path}

    timer = TimeReports(workspace=out_folder, options=run_options_m)

    result_init = RasterSplitter(wse_grid_folder=input_raster_grid_folder, output_folder=out_folder,
                                 tile_wse=tile_input_wse, split_by=split_raster_by,
                                 huc12_path=split_huc12_path, watch=timer)
    result_init.run_raster_splitter()

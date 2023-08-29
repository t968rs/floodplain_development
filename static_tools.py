import os
import arcpy


class StaticTools:

    @staticmethod
    def get_fema_sr():
        fema_sr = 'GEOGCS["GCS_North_American_1983",DATUM["D_North_American_1983",SPHEROID["GRS_1980",6378137.0,298.257222101]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],VERTCS["NAVD88_height_(ftUS)",VDATUM["North_American_Vertical_Datum_1988"],PARAMETER["Vertical_Shift",0.0],PARAMETER["Direction",1.0],UNIT["Foot_US",0.3048006096012192]];-400 -400 1000000000;-100000 10000;-100000 10000;8.98315284119521E-09;0.001;0.001;IsHighPrecision'

        return fema_sr

    @staticmethod
    def create_new_fema_dataset(output_gdb):

        arcpy.env.XYResolution = "0.0000000784415 DecimalDegrees"
        arcpy.env.XYTolerance = "0.0000000784415 DecimalDegrees"
        arcpy.env.ZTolerance = None
        if not arcpy.Exists(os.path.join(output_gdb, 'FIRM_Spatial_Layers')):
            arcpy.management.CreateFeatureDataset(out_dataset_path=output_gdb,
                                                  out_name='FIRM_Spatial_Layers',
                                                  spatial_reference=StaticTools.get_fema_sr())



class FEMAtables:

    def __init__(self, table_type):
        self.table_type = table_type

    @staticmethod
    def get_s_fld_field_infos():

        f_atts = ("Field", "R / A", "Type", "Length / Precision", "Joined", "Spatial / Lookup", "Domains")
        f_names = (
            'DFIRM_ID', 'VERSION_ID', 'FLD_AR_ID', 'STUDY_TYP', 'FLD_ZONE', 'ZONE_SUBTY', 'SFHA_TF', 'STATIC_BFE',
            'V_DATUM', 'DEPTH', 'LEN_UNIT', 'VELOCITY', 'VEL_UNIT', 'AR_REVERT', 'AR_SUBTRV', 'BFE_REVERT',
            'DEP_REVERT',
            'DUAL_ZONE', 'SOURCE_CIT')
        f_req = ('R', 'R', 'R', 'R', 'R', 'A', 'R', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'R')
        f_types = (
            'Text', 'Text', 'Text', 'Text', 'Text', 'Text', 'Text', 'Double', 'Text', 'Double', 'Text', 'Double',
            'Text',
            'Text', 'Text', 'Double', 'Double', 'Text', 'Text')
        length_prec = (6, 11, 25, 38, 17, 72, 1, '', 17, '', 16, '', 20, 17, 72, '', '', 1, 11)
        domain_relate = (
            '', '', '', 'D_Study_Typ', 'D_Zone', 'D_Zone_Subtype', 'D_TrueFalse', '', 'D_V_Datum', '', 'D_Length_Units',
            '',
            'D_Velocity_Units', 'D_Zone', 'D_Zone_Subtype', '', '', 'D_TrueFalse', 'L_Source_Cit')

        # Populate the list of dictionaries by iterating field names and the above tuple indices
        allfields = []
        for i, f_name in enumerate(f_names):
            thisfield_dict = {'field_name': f_name,
                              'field_type': f_types[i],
                              'Required': f_req[i]}
            if f_types[i] == 'Text':
                thisfield_dict['field_length'] = length_prec[i]
            else:
                thisfield_dict['field_length'] = ''
            if domain_relate[i] != '':
                if "D_" in domain_relate[i]:
                    thisfield_dict['field_domain'] = domain_relate[i]
                elif "L_" in domain_relate[i]:
                    thisfield_dict['Related Table'] = domain_relate[i]
                else:
                    thisfield_dict['field_domain'] = ''
                    thisfield_dict['Related Table'] = ''

            print(f'{f_name}\n  {thisfield_dict}')
            allfields.append(thisfield_dict)

        return tuple(allfields)

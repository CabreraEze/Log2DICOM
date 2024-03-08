"""DICOM related data structure"""
import pydicom

from .linac_class import AxisData, FluenceData

import numpy as np
import copy



class DicomBase:

    def __init__(self, file_path: str = None):
        self.raw_data = pydicom.dcmread(file_path)
        self.num_beams = len(self.raw_data.BeamSequence)
        self.control_points = self.get_cpoints()
        self.num_patient_beams = len(self.control_points)
        self.patientID = self.raw_data.PatientID

    def beam_order(self):
        beam_index = np.array([idx for idx in range(self.num_beams) if len(self.raw_data.BeamSequence[idx].ControlPointSequence) > 5])

        idx_order = []
        for beam in beam_index:
            order = 0
            while not 'field' + str(order + 1) in str(self.raw_data.BeamSequence[beam].BeamName).lower().replace(' ',''):
                order += 1
                # if order == 100: popAdvice('Invalid beam names'); exit()
            idx_order.append(order)
        return beam_index[idx_order].tolist()

    def get_cpoints(self):
        return [len(self.raw_data.BeamSequence[beam].ControlPointSequence) for beam in self.beam_order()]


# Main class for Dicom data
class DicomData(DicomBase):
    
    def __init__(self, file_path: str, qa_params: dict, equal_aspect: bool = True, is_hdmlc: bool = True):
        super().__init__(file_path)
        self.qa_params = qa_params
        self.axis = AxisData(
            self.raw_data,
            beams=self.beam_order(),
            control_points=self.control_points
        )
        self.fluence_map = FluenceData(self.axis, resolution=qa_params.get('resolution'), equal_aspect=equal_aspect, is_hdmlc=is_hdmlc)


class ModedDicom:

    def __init__(self, dicom, log):
        self.dicom_original = dicom.raw_data
        self.log_original = log.raw_data

        self.perform_changes(dicom, log)

    def perform_changes(self, dicom, log):
        plan = copy.deepcopy(self.dicom_original)
        beams = dicom.beam_order()

        prep_leaves = []
        leaves_data = copy.deepcopy(log.axis.leaf)
        for matrix in leaves_data:
            matrix = np.fliplr(matrix.transpose())
            matrix[:,:60], matrix[:,60:] = -np.fliplr(matrix[:,:60]), np.fliplr(matrix[:,60:])
            prep_leaves.append(matrix)

        for beam in beams:
            for point in range(1,dicom.control_points[beams.index(beam)]):
                for leaf in range(dicom.axis._num_leaves):
                    plan.BeamSequence[beam].ControlPointSequence[point].BeamLimitingDevicePositionSequence[0].LeafJawPositions[leaf] = np.round(prep_leaves[beams.index(beam)][point-1,leaf] * 10, decimals=1)
                    if hasattr(plan.BeamSequence[beams[0]].ControlPointSequence[1], 'GantryAngle'):
                        plan.BeamSequence[beam].ControlPointSequence[point].GantryAngle = log.axis.clockwise_counterclockwise(log.axis.gantry[beams.index(beam)][point-1], offset=180)
        
        for beam in range(len(plan.BeamSequence)):
            plan.BeamSequence[beam].BeamName += ' mod'
        plan.ApprovalStatus = 'UNAPPROVED'

        # Rename the UID, needed for import both the original and the modified plan to eclipse
        new_uid = plan.SOPInstanceUID.split('.')
        new_uid[-1] = str(int(new_uid[-1]) +1)
        new_uid = '.'.join(new_uid)
        plan.SOPInstanceUID = new_uid

        self.moded_plan = plan


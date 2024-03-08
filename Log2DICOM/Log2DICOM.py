"""Log2DICOM"""
import numpy as np
import matplotlib as mpl
import matplotlib.pylab as plt
import os
import argue

from core.dicom_class import DicomData, ModedDicom
from core.log_class import LogData
from core.tool import gamma_Bakai

class Log2DICOM:
    
    def __init__(
            self, log_paths, rt_path, qa_params,
            _is_hdmlc: bool = True, _fluence_w_points: bool = True, _equal_aspect: bool = True, _directory: str = os.getcwd()
        ):
        self._log_path = log_paths
        self._rt_path = rt_path
        self._qa_params = qa_params
        self._is_hdmlc = _is_hdmlc
        self._fluence_w_points = _fluence_w_points

        self.equal_aspect = _equal_aspect

        self._log_names = self.get_file_name('log')
        self._rt_name = self.get_file_name('rt')

        self.dicom = DicomData(self._rt_path, self._qa_params, equal_aspect=self.equal_aspect, is_hdmlc=self._is_hdmlc)
        self.log = LogData(self._log_path, self._qa_params, equal_aspect=self.equal_aspect, is_hdmlc=self._is_hdmlc, control_points=self.dicom.control_points, dicom_meterset=self.dicom.axis.mu)

        self.modified = ModedDicom(self.dicom, self.log)

        self.create_dir(_directory, self.dicom.patientID)

    @argue.options(file=('log', 'rt'))
    def get_file_name(self, file:str):
        file_path = getattr(self, '_' + file + '_path')
        return [fl.split('/')[-1] for fl in file_path]

    def public_results(self):


        def save_map(name: str, map):
            color = mpl.cm.Greys
            color = color.reversed()
            plt.imsave(fname=name, arr=map, cmap=color, format='png')


        def save_results(beam_name: str, log_fluence, dicom_fluence):
            save_map(self.patient_folder + '/log_' + beam_name + '_fluence.png', log_fluence)
            save_map(self.patient_folder + '/dicom_' + beam_name + '_fluence.png', dicom_fluence)

            gamma, gamma_report = gamma_Bakai(log_fluence, dicom_fluence, report=True, **self._qa_params)
            save_map(self.patient_folder + '/gamma_' + beam_name + '.png', gamma)

            with open(self.patient_folder + '/gamma_'+ beam_name +'_report.txt', 'w') as f:
                print(gamma_report, file=f)


        if not self._fluence_w_points and len(self._log_names)!=self.dicom.num_beams:
            log_fluence = self.log.raw_data[0].fluence.actual.calc_map(equal_aspect=self.equal_aspect)
            dicom_fluence = np.array(self.dicom.fluence_map.array).sum(axis=0)
            save_results('total', log_fluence, dicom_fluence)
        else:
            if not self._fluence_w_points:
                log_fluence = [self.log.raw_data[beam].fluence.actual.calc_map(equal_aspect=self.equal_aspect) for beam in range(self.log.axis._num_beams)]
            else:
                log_fluence = self.log.fluence_map.array

            for beam in range(self.log.axis._num_beams):
                dicom_fluence = self.dicom.fluence_map.array[beam]
                save_results('field'+str(beam), log_fluence[beam], dicom_fluence)

    def save(self):
        self.modified.moded_plan.save_as(self.patient_folder + '/' + self.dicom.patientID + '_mod.dcm')

    def create_dir(self, directory:str, name:str):
        self.patient_folder = directory + '/' + name
        if not os.path.exists(self.patient_folder):
            os.makedirs(self.patient_folder)
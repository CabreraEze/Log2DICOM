"""LogFile related data structure"""
from pylinac.log_analyzer import load_log

from .tool import snapshot_to_point
from .linac_class import AxisData, FluenceData

from datetime import datetime


class LogBase:

    def __init__(self, log_path: list):
        self._num_logs = len(log_path)
        self._log_path = log_path
        self.raw_data = self._load_log()
    
    def _load_log(self):
        log_list = list(self._log_path)

        def log_by_field(path):
            order = 0
            while not 'field' + str(order+1) in path.split('/')[-1].lower().replace(' ',''):
                order += 1
                # if order == 100: popAdvice('Invalid logfiles names'); exit()
            return order
        

        if log_list[0].split('.')[-1] == 'bin':
            log_list.sort(key=log_by_field)
        else:
            log_list.sort(key=lambda path: datetime.strptime(path.split('/')[-1].split('.')[0].split('_')[0][1:], '%Y%m%d%H%M%S'))
        return [load_log(path) for path in log_list]


# Main class for Logs data
class LogData(LogBase):
    
    def __init__(self, log_path: list, qa_params: dict, equal_aspect: bool = True, is_hdmlc: bool = True, control_points: list = None, dicom_meterset: list = None):
        super().__init__(log_path)
        self._dicom_points = control_points
        self._snap_for_point = snapshot_to_point(logs=self.raw_data, control_points=control_points, meterset=dicom_meterset)
        self.axis = AxisData(self.raw_data, nlogs=self._num_logs, control_points=self._dicom_points, snapshots=self._snap_for_point)
        self.fluence_map = FluenceData(self.axis, resolution=qa_params.get('resolution'), equal_aspect=equal_aspect, is_hdmlc=is_hdmlc)


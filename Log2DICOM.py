#############################################################################
#                                                                           #
#   Log2DICOM: A python-based tool for Modifying original DICOM planning    #
#   files based on measured treatment data                                  #
#                                                                           #
#   This script extracts relevant information from log files and uses that  #
#   data to make appropriate modifications to DICOM files. Designed with    #
#   precision and user-friendliness in mind, it streamlines the process of  #
#   adjusting medical imaging data based on logged operations or events.    # 
#                                                                           #
#   Author: Ezequiel Agustin Cabrera Lacey                                  #
#   Date: 9-13-2023                                                         #
#   Version: 0.1                                                            #      
#                                                                           #
#   Recognitions:                                                           #
#   -National University of Córdoba an Zunino Institute to provide the      #
#    knowledge and space to make this project posible.                      #  
#   -Dr. Miguel Chesta and Dr. Daniel Venencia for their contributions,     #
#    insights, and support during the development of this scritp.           #
#                                                                           #
#   DISCLAIMER:                                                             #
#   This software constitude the final project of a physics degree and is   #
#   provided "as is", without warranty of any kind, express or implied,     #
#   including but not limited to the warranties of merchantability, fitness #
#   for a particular purpose, and non-infringement.                         #
#   In no event shall the authors or copyright holders be liable for any    #
#   claim, damages or other liability, whether in an action of contract,    #
#   tort or otherwise, arising from, out of or in connection with the       #
#   software or the use or other dealings in the software. Users are        #
#   encouraged to test the software extensively before relying on it for    #
#   critical applications.                                                  #
#                                                                           #
#############################################################################




from pylinac.log_analyzer import load_log
from pylinac import TrajectoryLog, Dynalog, MachineLogs
from pydicom.dataset import FileDataset
from scipy.ndimage import sobel
import pydicom

from tkinter import filedialog, Button, Label, Entry
from sys import exit
from bisect import bisect_left, bisect_right

import tkinter as tk
import numpy as np
import matplotlib as mpl
import matplotlib.pylab as plt
import copy
import os
import argue




class mainWindow():
    """
    Main class for the UI interface, mainly for loading DICOM and Log files path and setting
    gamma analysis parameters

    Atributes:
        distTA (float): gamma analysis distance to agreatment
        doseTA (float): gamma analysis dose to agreatment
        threshold (float): gamma analysis threshold for calculation
        resolution (float): fluence map resolution 

        -not accesible-
        _is_hd_mlc (bool): whether or not, the mlc has high definition
        _fluence_w_points (bool): whether or not the log fluence map is made only with the control point data
        _directory (str): directory where the modified DICOM file will be located
        _equal_aspect (bool): matplotlib equal aspect parameter
    """

    def __init__(self,master):
        self.master = master

        self.rt_button = Button(master, text='Load Dicom plan', command=self.load_dicom)
        self.rt_button.pack()

        self.log_button = Button(master, text='Load Log\s', command=self.load_log)
        self.log_button.pack()

        self.next_button = Button(master, text='Next', command=self.QA_Params, background='orange')
        self.next_button.pack()

    def load_log(self):
        self.log_paths = filedialog.askopenfilenames(title='Import LogFile\s', filetypes=(('Dynalog files', '*.dlg'), ('Trajectorylog files', '*.bin'), ('all files', '*.*')))

    def load_dicom(self):
        self.rt_path = filedialog.askopenfilename(title='Import Dicom plan', filetypes=(('Dicom files', '*.dcm'), ('all files', '*.*')))

    def QA_Params(self):

        self.rt_button.destroy()
        self.log_button.destroy()
        self.next_button.destroy()

        # default qa params
        ###########################
        self.distTA = 1
        self.doseTA = 1
        self.resolution = 0.1
        self.threshold = 0.1
        ###########################

        def validate(text:str):
            return text.replace('.','',1).isdecimal() or text==''
        
        self.distTA_label = Label(self.master, text='Distance TA [mm]:')
        self.distTA_entry = Entry(
            self.master,
            validate='key',
            validatecommand=(self.master.register(validate), '%P')
        )
        self.distTA_label.pack()
        self.distTA_entry.pack()

        self.doseTA_label = Label(self.master, text='Dose TA [%]:')
        self.doseTA_entry = Entry(
            self.master,
            validate='key',
            validatecommand=(self.master.register(validate), '%P')
        )
        self.doseTA_label.pack()
        self.doseTA_entry.pack()

        self.resolution_label = Label(self.master, text='Fluence resolution [mm]:')
        self.resolution_entry = Entry(
            self.master,
            validate='key',
            validatecommand=(self.master.register(validate), '%P')
        )
        self.resolution_label.pack()
        self.resolution_entry.pack()

        self.threshold_label = Label(self.master, text='Gamma threshold [%]:')
        self.threshold_entry = Entry(
            self.master,
            validate='key',
            validatecommand=(self.master.register(validate), '%P')
        )
        self.threshold_label.pack()
        self.threshold_entry.pack()

        ###########################
        # Change to be customizable
        self._is_hdmlc = True
        self._fluence_w_points = True
        self._directory = os.getcwd()
        self._equal_aspect = True
        ###########################

        self.default_val = False
        self.default_val_buttom = Button(self.master, text='Set standard values', command=self.get_def_val, background='cyan')
        self.default_val_buttom.pack()

        self.entry = Button(self.master, text='Next', command=self.get_values, background='orange')
        self.entry.pack()
    
    def get_def_val(self):
        self.default_val = True
        self.get_values()

    def get_values(self):
        if not self.default_val:
            self.distTA = float(self.distTA_entry.get())
            self.doseTA = float(self.doseTA_entry.get())
            self.resolution = float(self.resolution_entry.get())
            self.threshold = float(self.threshold_entry.get())/100
        self.master.destroy()

    @property
    def get_params(self):
        list_params = [
            'distTA',
            'doseTA',
            'resolution',
            'threshold'
        ]
        params = {}
        for param in list_params:
            params[param] = getattr(self,param)
        return params




class popAdvice:
    """
    Poping warnings when needing:

    Types of warnings:
        'Invalid beam names': Dicom field names must have 'Field x' in their name where x=field number
        'Invalid logfiles names': For multiple field, Dynalog files must have field specified in their name
        '¨kind¨ doesnt math any aceptable case': For interpolate TrajectoryLog file's control point snapshot to DICOM's control points,
            must select 'low', 'high' or 'medium'
    """

    def __init__(self ,advice: str):
        self.master = tk.Tk()

        self.advice = Label(self.master, text=advice)
        self.advice.pack()

        self.exit_buttom = Button(self.master, text='Close', command=self.master.destroy)
        self.exit_buttom.pack()

        self.master.mainloop()




class Log2DICOM:
    """
    Main class for performing DICOM modifications, for its use import DICOM rt plan and its corresponding
    Varian Logfile/s. This class perform gamma analysis on the reconstructed fluence.
    The results for the plan modifications and gamma analysis will be exported to a folder with the patient ID as its name

    Atributes:
        _log_path (list): List of directions where the Log files would be located
        _rt_path (str): String containing the direction where the DICOM plan would be loaded
        _qa_params (dict): Dictionary with the gamma analysis parameters
        _is_hdmlc (bool): Whether the mlc has high definition
        _fluence_w_points (bool): whether or not the log fluence map is made only with the control point data
        equal_aspect (bool): matplotlib equal aspect parameter
        _log_names (list): Log file/s name/s
        _rt_name (str): DICOM file name

        dicom (DicomData): Class containing the recovered data from the rt plan
        log (LogData): Class containing the recovered data from the corresponding log/s
        modified (ModedDicom): Class containing the modified DICOM plan

    Methods:
        public_results: Save the results from the gamma analysis
        save: Save the modified DICOM plan
    """

    def __init__(self, user_loaded: mainWindow):
        self._log_path = user_loaded.log_paths
        self._rt_path = user_loaded.rt_path
        self._qa_params = user_loaded.get_params
        self._is_hdmc = user_loaded._is_hdmlc
        self._fluence_w_points = user_loaded._fluence_w_points

        self.equal_aspect = user_loaded._equal_aspect

        self._log_names = self.get_file_name('log')
        self._rt_name = self.get_file_name('rt')

        self.dicom = DicomData(self._rt_path, self._qa_params, equal_aspect=self.equal_aspect, is_hdmlc=self._is_hdmc)
        self.log = LogData(self._log_path, self._qa_params, equal_aspect=self.equal_aspect, is_hdmlc=self._is_hdmc, control_points=self.dicom.control_points, dicom_meterset=self.dicom.axis.mu)

        self.modified = ModedDicom(self.dicom, self.log)

        self.create_dir(user_loaded._directory, self.dicom.patientID)

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




class DicomBase:
    """
    Base class with basic DICOM information

    Atributes:
        raw_data (pydicom.FileDataSet): Pydicom DICOM object
        num_beams (int): Number of beams (could not be corresponded with the number of treatment Fields)
        control_points (list): List containing the number of control points in each beam
        num_patient_beams (int): Number of Fields (could not be corresponded with the number of beams)
        patiendID (str): PatiendID

    Methods:
        beam_order: Return the beam number corresponding with each field, ordered by name
    """

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
            while 'Field ' + str(order +1) != self.raw_data.BeamSequence[beam].BeamName:
                order += 1
                if order == 100: popAdvice('Invalid beam names'); exit()
            idx_order.append(order)
        return beam_index[idx_order].tolist()

    def get_cpoints(self):
        return [len(self.raw_data.BeamSequence[beam].ControlPointSequence) for beam in self.beam_order()]




class DicomData(DicomBase):
    """
    Main class for storing DICOM data

    Atributes:
        qa_params (dict): Dictionary with the gamma analysis parameters
        axis (AxisData): Stored leaf, jaw, gantry and mu (meterset for DICOM) data
        fluence_map (FluenceData): Reconstructed fluence map
    """

    def __init__(self, file_path: str, qa_params: dict, equal_aspect: bool = True, is_hdmlc: bool = True):
        super().__init__(file_path)
        self.qa_params = qa_params
        self.axis = AxisData(
            self.raw_data,
            beams=self.beam_order(),
            control_points=self.control_points
        )
        self.fluence_map = FluenceData(self.axis, resolution=qa_params.get('resolution'), equal_aspect=equal_aspect, is_hdmlc=is_hdmlc)




class AxisData:
    """
    Main class for storing the axis control point data in Log file format 

    Atributes:
        _from_type (type): Type of data where the axis is retrieved
        ganry (list): Gantry angle in each control point for each beam
        jaws (JawData): Stored leaf, jaw, gantry and mu (meterset for DICOM) data
        leaf (list): List with the position of each leaf on each snapshot, for each treatment beam
        mu (list): Stored; Monitor Units for TrajectoryLog, Normalized dose fraction for Dynalog and Meterset for Dicom files, for each beam

    Methods:
        clockwise_counterclockwise: Perfom radial coordinates transformations
    """

    def __init__(self, data_object, **kwargs):   
        self._from_type = type(data_object)

        self.gantry = self.get_axis(data_object, 'gantry', **kwargs)
        self.jaws = self.get_axis(data_object, 'jaws', **kwargs)
        self.leaf = self.get_axis(data_object, 'leaf', **kwargs)
        self.mu = self.get_axis(data_object, 'mu', **kwargs)

    axis_options = ('leaf', 'gantry', 'jaws', 'mu')
    def get_axis(self, data, axis: str, **kwargs):
        if isinstance(data, FileDataset):
            return self.get_from_dicom(data, axis, **kwargs)
        else:
            return self.get_from_log(data, axis, **kwargs)

    @argue.options(axis=axis_options)
    def get_from_dicom(self, data, axis: str, beams: list = None, control_points: list = None):
        self._num_beams = len(beams)
        self._control_points = control_points
        if axis == 'jaws':
            return JawData(data, beams=beams)
        else:
            return self.axis_from_dicom(data, axis, beams, control_points)

    @argue.options(axis=axis_options)
    def get_from_log(self, data, axis: str, nlogs: int = None, control_points: list = None, snapshots: list = None):
        self._num_beams = len(control_points)
        self._control_points = control_points
        if axis == 'jaws':
            return JawData(data, nlogs=nlogs, control_points=control_points , snapshots=snapshots)
        else:
            return self.axis_from_log(data, axis, control_points=control_points, snapshots=snapshots)

    @argue.options(axis=axis_options)
    def axis_from_dicom(self, dicom: FileDataset, axis: str, beams: list, control_points: list):
        if axis == 'leaf':
            leaf_axes = []
            for beam in beams:
                leaf_data = np.transpose(
                            [dicom.BeamSequence[beam].ControlPointSequence[point].BeamLimitingDevicePositionSequence[0].LeafJawPositions \
                            for point in range(1,control_points[beams.index(beam)])],
                            axes=(1,0)
                )
                
                # For consistency with the logs
                leaf_data = np.flipud(leaf_data)
                leaf_data[:60,:], leaf_data[60:,:] = np.flipud(leaf_data[:60,:]), -np.flipud(leaf_data[60:,:])
                leaf_axes.append(leaf_data / 10)
            self._num_leaves = np.shape(leaf_axes[0])[0]
            return leaf_axes
        elif axis == 'gantry':
            gantry = []
            for beam in beams:
                gantry_data = [self.clockwise_counterclockwise(dicom.BeamSequence[beam].ControlPointSequence[point].GantryAngle, offset=180) for point in range(1, control_points[beams.index(beam)])]
                gantry.append(gantry_data)
            return gantry
        else:
            mu = []     
            for beam in beams: 
                mu_data = [dicom.BeamSequence[beam].ControlPointSequence[point].ReferencedDoseReferenceSequence[0].CumulativeDoseReferenceCoefficient for point in range(1, control_points[beams.index(beam)])]
                mu.append(mu_data)
            return mu

    @argue.options(axis=axis_options)
    def axis_from_log(self, data: list, axis: str, control_points: list = None, snapshots: list = None):
        num_beams = len(control_points)
        self._num_leaves = data[0].header.num_mlc_leaves
        nlogs = len(data)

        @argue.options(axis=('leaf', 'gantry', 'mu'))
        def _get_axis(axis: str, log, snapshots: list):
            snapshots = snapshots[1:]
            if axis == 'leaf':
                log_pos = [log.axis_data.mlc.leaf_axes[1].actual[snap] for snap in snapshots]
                for leaf in range(1, self._num_leaves):
                    log_pos = np.vstack(
                        (log_pos, [log.axis_data.mlc.leaf_axes[leaf+1].actual[snap] for snap in snapshots])
                    )
                return log_pos
            elif axis == 'gantry':
                return [log.axis_data.gantry.actual[snap] for snap in snapshots]
            else:
                return [log.axis_data.mu.actual[snap] for snap in snapshots]

        if nlogs != num_beams:
            return[_get_axis(axis=axis, log=data[0], snapshots=snapshots[beam]) for beam in range(num_beams)]
        else:
            return [_get_axis(axis=axis, log=log, snapshots=snapshots[beam]) for log, beam in zip(data, range(num_beams))]

    # For consistency with the logs
    def clockwise_counterclockwise(self, theta, offset: float = 0):
        return (-theta + offset) % 360




class FluenceData:
    """
    Class for storing reconstructed fluence map

    Atributes:
        array (np.array): Map array
    """

    def __init__(self, axis: AxisData = None, resolution: float = 0.1, equal_aspect: bool = True, is_hdmlc: bool = True):
        self.is_hdmlc = is_hdmlc
        self.resolution = resolution
        self.equal_aspect = equal_aspect
        self.axis = axis

        self.array: np.array = self.calc_map()

    def calc_map(self):
        resolution = self.resolution
        is_hdmlc = self.is_hdmlc
        equal_aspect = self.equal_aspect

        MLC_FOV_HEIGHT_MM = 400
        HDMLC_FOV_HEIGHT_MM = 220
        MLC_FOV_WIDTH_MM = 400

        num_beams = self.axis._num_beams
        control_points = self.axis._control_points
        num_pairs = self.axis._num_leaves // 2
        leaf_axes = self.axis.leaf

        x1 = self.axis.jaws.x1
        y1 = self.axis.jaws.y1

        x2 = self.axis.jaws.x2
        y2 = self.axis.jaws.y2

        height = MLC_FOV_HEIGHT_MM if not is_hdmlc else HDMLC_FOV_HEIGHT_MM
        if equal_aspect:
            fluence = np.zeros(
                (num_beams, int(height / resolution), int(MLC_FOV_WIDTH_MM / resolution)),
                dtype=float,
            )
        else:
            fluence = np.zeros(
                (num_beams, num_pairs, int(MLC_FOV_WIDTH_MM / resolution)), dtype=float
            )

        def create_mlc_y_positions(is_hdmlc):
            if not is_hdmlc:
                num_large_leaves = 10
                size_large_leaves = 10 / resolution
                num_small_leaves = 40
                size_small_leaves = 5 / resolution
            else:
                num_large_leaves = 14
                size_large_leaves = 5 / resolution
                num_small_leaves = 32
                size_small_leaves = 2.5 / resolution
            sizes = (
                [size_large_leaves] * num_large_leaves
                + [size_small_leaves] * num_small_leaves
                + [size_large_leaves] * num_large_leaves
            )
            return np.cumsum(
                [
                    0,
                ]
                + sizes
            ).astype(int)

        positions = create_mlc_y_positions(is_hdmlc=is_hdmlc)

        def yield_leaf_width():
            for idx in range(num_pairs):
                yield (positions[idx], positions[idx + 1])

        # calculate the MU delivered in each snapshot. For Tlogs this is absolute; for dynalogs it's normalized.
        mu_matrix = self.axis.mu
        MU_differential = [np.array([mu_matrix[beam][0]] + list(np.diff(mu_matrix[beam][:]))) for beam in range(num_beams)]
        MU_total = [mu_matrix[beam][-1] for beam in range(num_beams)]

        # calculate each "line" of fluence (the fluence of an MLC leaf pair, e.g. 1 & 61, 2 & 62, etc),
        # and add each "line" to the total fluence matrix
        def leaf_under_y_jaw(
                self, leaf_num: int, y1: float, y2: float
            ) -> bool:
            """Return a boolean specifying if the given leaf is under one of the y jaws.

            Parameters
            ----------
            is_hdmlc: bool
            leaf_num: int
            y1 : float
            y2 : float
            """
            # Changes in base of varian mlc models
            inner_leaf = 40
            outer_leaf = 10

            outer_leaf_thickness = 10  # mm
            inner_leaf_thickness = 5
            mlc_position = 0
            if self.is_hdmlc:
                # Varian hd mlc models have less outer leaves
                inner_leaf = 32
                outer_leaf = 14

                outer_leaf_thickness /= 2
                inner_leaf_thickness /= 2
                mlc_position = 90

            # Changes performed to reflect the differences between hd and standard varian mlc models
            for leaf in range(1, leaf_num + 1):
                if outer_leaf >= leaf or leaf >= 120 - outer_leaf:
                    mlc_position += outer_leaf_thickness
                elif inner_leaf + outer_leaf >= leaf or leaf >= 60 + outer_leaf:
                    mlc_position += inner_leaf_thickness
                else:  # in between
                    mlc_position += outer_leaf_thickness

            y2_position = y2 * 10 + 200
            y1_position = 200 - y1 * 10
            if 10 >= leaf_num or leaf_num >= 110:
                thickness = outer_leaf_thickness
            elif 50 >= leaf_num or leaf_num >= 70:
                thickness = inner_leaf_thickness
            else:  # in between
                thickness = outer_leaf_thickness
            return mlc_position < y1_position or mlc_position - thickness > y2_position

        def pair_moved(pair_num: int, beam:int) -> bool:
            """Return whether the given pair moved during treatment.

            If either leaf moved, the pair counts as moving.

            Parameters
            ----------
            pair_num : int
            beam : int
            """
            a_leaf = pair_num
            b_leaf = pair_num + num_pairs
            return leaf_moved(a_leaf, beam) or leaf_moved(b_leaf, beam)

        def leaf_moved(leaf_num: int, beam:int) -> bool:
            """Return whether the given leaf moved during treatment.

            Parameters
            ----------
            leaf_num : int
            beam : int
            """
            return leaf_num in moving_leaves(beam)

        def moving_leaves(beam:int) -> np.ndarray:
            """Return an array of the leaves that moved during treatment."""
            threshold = 0.01
            indices = []
            for leaf_num, leafdata in zip(range(num_pairs*2), leaf_axes[beam]):
                leaf_stdev = np.std(leafdata)
                if leaf_stdev > threshold:
                    indices.append(leaf_num)
            return indices


        fluence_line = [np.zeros((int(400 / resolution)), dtype=np.float32) for beam in range(num_beams)]
        pos_offset = int(np.round(200 / resolution))
        for beam in range(num_beams):
            for pair, width in zip(range(num_pairs), yield_leaf_width()):
                if not leaf_under_y_jaw(pair, y1=y1[beam], y2=y2[beam]):
                    fluence_line[beam][:] = 0  # emtpy the line values on each new leaf pair

                    right_leaf_data = leaf_axes[beam][pair, :]
                    right_leaf_data = (
                        np.round(right_leaf_data * 10 / resolution) + pos_offset
                    )
                    left_leaf_data = leaf_axes[beam][pair + num_pairs, :]

                    left_leaf_data = (
                        -np.round(left_leaf_data * 10 / resolution) + pos_offset
                    )
                    left_jaw_data = np.round(
                        (200 / resolution) - (x1[beam] * 10 / resolution)
                    )
                    right_jaw_data = np.round(
                        (x2[beam] * 10 / resolution) + (200 / resolution)
                    )
                    if pair_moved(pair, beam):
                        for snapshot in range(control_points[beam]-1):
                            lt_mlc_pos = left_leaf_data[snapshot]
                            rt_mlc_pos = right_leaf_data[snapshot]
                            lt_jaw_pos = left_jaw_data
                            rt_jaw_pos = right_jaw_data
                            left_edge = int(max(lt_mlc_pos, lt_jaw_pos))
                            right_edge = int(min(rt_mlc_pos, rt_jaw_pos))
                            fluence_line[beam][left_edge:right_edge] += MU_differential[beam][snapshot]
                    else:  # leaf didn't move; no need to calc over every snapshot
                        first_snapshot = 0
                        lt_mlc_pos = left_leaf_data[first_snapshot]
                        rt_mlc_pos = right_leaf_data[first_snapshot]
                        lt_jaw_pos = left_jaw_data
                        rt_jaw_pos = right_jaw_data
                        left_edge = max(lt_mlc_pos, lt_jaw_pos)
                        right_edge = min(rt_mlc_pos, rt_jaw_pos)
                        fluence_line[beam][int(left_edge) : int(right_edge)] = MU_total[beam]
                    if equal_aspect:
                        fluence[beam][width[0] : width[1], :] = np.tile(
                            fluence_line[beam], [width[1] - width[0], 1]
                        )
                    else:
                        fluence[beam][pair - 1, :] = fluence_line[beam]

        #fluence[beam][y_axis][x_axis]
        #planned_fluence[beam]
        self.planned_fluence = MU_total
        return fluence




class JawData:
    """
    Stored jaw data in Log format

    Atributes
        x1 (list): X1 jaw position for each field
        x2 (list): X2 jaw position for each field
        y1 (list): Y1 jaw position for each field
        y2 (list): Y2 jaw position for each field
    """
        
    def __init__(self, data, **kwargs):

        x_jaw, y_jaw = self.jaw_from_data(data, **kwargs)

        self.x1 = x_jaw[:,0]
        self.y1 = y_jaw[:,0]

        self.x2 = x_jaw[:,1]
        self.y2 = y_jaw[:,1]

    def jaw_from_dicom(self, data, beams):
        x_jaw = [
            np.array(data.BeamSequence[beam].ControlPointSequence[0].BeamLimitingDevicePositionSequence[0].LeafJawPositions) \
            for beam in beams
        ]
        y_jaw = [
            np.array(data.BeamSequence[beam].ControlPointSequence[0].BeamLimitingDevicePositionSequence[1].LeafJawPositions) \
            for beam in beams
        ]
        # For consistency with the logs
        x_jaw = np.array([jaw/10 for jaw in x_jaw])
        y_jaw = np.array([jaw/10 for jaw in y_jaw])

        x_jaw[:,0] = -x_jaw[:,0]
        y_jaw[:,0] = -y_jaw[:,0]
        return x_jaw, y_jaw

    def jaw_from_data(self, data, **kwargs):
        if isinstance(data, FileDataset):
            return self.jaw_from_dicom(data, beams=kwargs.get('beams'))
        else:
            return self.jaw_from_log(data, **kwargs)

    def jaw_from_log(self, logs, nlogs: int = None, control_points: list = None, snapshots: list = None):
        num_beams = len(control_points)
        if nlogs != num_beams:
            x_jaw = [[
                logs[0].axis_data.jaws.x1.actual[snapshots[beam][1]],
                logs[0].axis_data.jaws.x2.actual[snapshots[beam][1]]
            ] for beam in range(num_beams)] 

            y_jaw = [[
                logs[0].axis_data.jaws.y1.actual[snapshots[beam][1]],
                logs[0].axis_data.jaws.y2.actual[snapshots[beam][1]]
             ] for beam in range(num_beams)]      
        else:
            x_jaw = []
            y_jaw = []
            for log in logs:
                x_jaw.append([
                    log.axis_data.jaws.x1.actual[0],
                    log.axis_data.jaws.x2.actual[0]
                ])
                
                y_jaw.append([
                    log.axis_data.jaws.y1.actual[0],
                    log.axis_data.jaws.y2.actual[0]
                ])
        return np.array(x_jaw), np.array(y_jaw)




class LogBase:
    """
    Base class for store Log basic data

    Atributes:
        _num_logs (int): Num of logs loaded
        _log_path (list): List containing the directory for each log loaded
        raw_data (pylinac.TrajectoryLog, pylinac.Dynalog): Pylinac Log object
    """

    def __init__(self, log_path: list):
        self._num_logs = len(log_path)
        self._log_path = log_path
        self.raw_data = self._load_log()
    
    def _load_log(self):
        log_list = list(self._log_path)

        def log_by_field(path):
            order = 0
            while not 'Field ' + str(order+1) in path:
                order += 1
                if order == 100: popAdvice('Invalid logfiles names'); exit()
            return order
        log_list.sort(key=log_by_field)
        return [load_log(path) for path in log_list]




class LogData(LogBase):
    """
    Main class for storing Log data

    Atributes:
        _dicom_points (list): List with the number of control points for each field
        _snap_for_points (list): List of snapshots indices to match each control point
        axis (AxisData): Stored axis data for each control point
        fluence_map (FluenceData): Reconstructed fluence map
    """

    def __init__(self, log_path: list, qa_params: dict, equal_aspect: bool = True, is_hdmlc: bool = True, control_points: list = None, dicom_meterset: list = None):
        super().__init__(log_path)
        self._dicom_points = control_points
        self._snap_for_point = snapshot_to_point(logs=self.raw_data, control_points=control_points, meterset=dicom_meterset)
        self.axis = AxisData(self.raw_data, nlogs=self._num_logs, control_points=self._dicom_points, snapshots=self._snap_for_point)
        self.fluence_map = FluenceData(self.axis, resolution=qa_params.get('resolution'), equal_aspect=equal_aspect, is_hdmlc=is_hdmlc)




class ModedDicom:
    """
    Main class for storing and perfom modifications in the original DICOM rt plan

    Atributes:
        dicom_original (pydicom.FileDataSet): Original DICOM rt plan
        log_original (pylinac.TrajectoryLog, pylinac.Dynalog): Log file data
        moded_plan (pydicom.FileDataSet): Stored DICOM modified plan
    """

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




def gamma_Bakai(
            reference: np.ndarray,
            comparison: np.ndarray,
            resolution: float = 0.1,
            threshold: float = 0.1,
            distTA: float = 1,
            doseTA: float = 1,
            report: bool = False
    ) -> np.ndarray:
    """
    Calculate the gamma from the actual and expected fluences.

    The gamma calculation is based on `Bakai et al
    <http://iopscience.iop.org/0031-9155/48/21/006/>`_ eq.6,
    which is a quicker alternative to the standard Low gamma equation.

    Parameters
    ----------
    reference: numpy.ndarray
        Reference map
    comparison: numpy.ndarray
        Comparison map
    doseTA : int, float
        Dose-to-agreement in percent; e.g. 2 is 2%.
    distTA : int, float
        Distance-to-agreement in mm.
    threshold : int, float
        The dose threshold percentage of the maximum dose, below which is not analyzed.
    resolution : int, float
        The resolution in mm of the resulting gamma map in the leaf-movement direction.
    calc_individual_maps : bool
        Not yet implemented.
        If True, separate pixel maps for the distance-to-agreement and dose-to-agreement are created.

    Returns
    -------
    numpy.ndarray
        A num_mlc_leaves-x-400/resolution numpy array.
    """

    ref_img = copy.deepcopy(reference)
    comp_img = copy.deepcopy(comparison)

    ref_img /= ref_img.max()
    comp_img /= comp_img.max()

    dpi = 25.4 / resolution
    dpmm = dpi / 25.4

    # invalidate dose values below threshold so gamma doesn't calculate over it
    ref_img[ref_img < threshold * np.max(ref_img)] = np.NaN

    # convert distance value from mm to pixels
    distTA_pixels = dpmm * distTA

    # construct image gradient using sobel filter
    img_x = sobel(ref_img.astype(np.float32), 1)
    img_y = sobel(ref_img.astype(np.float32), 0)
    grad_img = np.hypot(img_x, img_y)

    # equation: (measurement - reference) / sqrt ( doseTA^2 + distTA^2 * image_gradient^2 )
    subtracted_img = np.abs(comp_img - ref_img)
    denominator = np.sqrt(
        ((doseTA / 100.0) ** 2) + ((distTA_pixels**2) * (grad_img**2))
    )
    gamma_map = subtracted_img / denominator

    def GammaReport(gamma_map: np.ndarray):
        pixels_passing = np.sum(gamma_map[~np.isnan(gamma_map)] < 1)
        all_calcd_pixels = np.sum(gamma_map[~np.isnan(gamma_map)] >= 0)
        report = 'pixels_passing:\t\t' + str(pixels_passing) +\
        '\nall_calcd_pixels:\t' + str(all_calcd_pixels) +\
        '\npass_prcnt:\t\t' + str(np.round(pixels_passing / all_calcd_pixels * 100, 2)) +\
        '\nDTA:\t\t\t' + str(distTA) +\
        '\nDD:\t\t\t' + str(doseTA) +\
        '\nthreshold:\t\t' + str(threshold * 100)
        return report

    if report:
        gamma_report = GammaReport(gamma_map)
        return np.nan_to_num(gamma_map), gamma_report
    else:
        return np.nan_to_num(gamma_map)


def snapshot_to_point(
    logs, control_points: list, kind: str = 'low', meterset: list = None
    ) -> list:
    '''
    Set a correlation between Logfile snapshot and plan Control Points
    this only can be performed on TrayectoryLog files.
    From convenience, the snapshot-control point index will be consider
    the last measurement in that control point

    Parameters
    ----------
    logs: pylinac-object
        Axis object containing the control points
    kind: 'low', 'high', 'medium'
        The type of criterion to be consider the representative snapshot
        for each control point

    Return
    ---------
    indices: list
    '''


    def interpol_snap_dynalog(log, control_points: int, meterset: list) -> list:
        dyna_meterset = log.axis_data.mu.actual / 25000
        
        indices = [0,]
        for point_meterset in meterset:
            diff = dyna_meterset - point_meterset
            indices.append(bisect_left(diff,0))
        return indices


    @argue.options(kind=('low', 'high', 'medium'))
    def snap_knd(kind: str, control_points: int, log_snap: list, offset: int = 0):
        last_point = control_points + offset
        num_snap = len(log_snap)

        indices = []
        prev_point = 0
        if kind == 'low':
            for point in range(offset, last_point):
                for idx in range(prev_point, num_snap):
                    if log_snap[idx] == point:
                        prev_point = idx
                        indices.append(idx)
                        break
        elif kind == 'high':
            for point in range(offset, last_point):
                for idx in range(num_snap-1, prev_point-1, -1):
                    if log_snap[idx] == point:
                        prev_point = idx
                        indices.append(idx)
                        break
        elif kind == 'medium':
            low_idx = []
            for point in range(offset, last_point):
                for idx in range(prev_point, num_snap):
                    if log_snap[idx] == point:
                        prev_point = idx
                        low_idx.append(idx)
                        break
            
            high_idx = []
            for point in range(offset, last_point):
                for idx in range(num_snap-1, prev_point-1, -1):
                    if log_snap[idx] == point:
                        prev_point = idx
                        high_idx.append(idx)
                        break
            
            indices = [int(low_idx + 0.5*(high_idx[snap] - low_idx[snap])) for snap in range(len(low_idx))]
        else:
            raise popAdvice('¨kind¨ doesnt math any aceptable case')
        return indices


    num_beams = len(control_points)
    nlogs = len(logs)
    if isinstance(logs[0], (TrajectoryLog,)):
        log_control_points = [[int(point) for point in log.axis_data.control_point.actual] for log in logs]
    else:
        return [interpol_snap_dynalog(logs[beam], control_points[beam], meterset[beam]) for beam in range(nlogs)]

    if nlogs != num_beams:
        offset = 0
        indices = []
        for beam in range(num_beams):
            indices.append(snap_knd(
                kind=kind,
                control_points=control_points[beam],
                log_snap=log_control_points[0],
                offset=offset
            ))
            offset += control_points[beam]
        return indices
    else:
        return [snap_knd(
            kind=kind,
            control_points=control_points[beam],
            log_snap=log_control_points[beam]
        ) for beam in range(num_beams)]
    



if True:
    root = tk.Tk()
    root.geometry('320x240')
    root.title('Comparative Wizard_demo')
    main = mainWindow(root)
    root.mainloop()
    
    sesion = Log2DICOM(main)
    print('saving...')
    sesion.save()
    sesion.public_results()
"""Linac related data structure"""
from pydicom.dataset import FileDataset

import numpy as np
import argue


class AxisData:

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
            if hasattr(dicom.BeamSequence[beams[0]].ControlPointSequence[1], 'GantryAngle'):
                gantry = []
                for beam in beams:
                    gantry_data = [self.clockwise_counterclockwise(dicom.BeamSequence[beam].ControlPointSequence[point].GantryAngle, offset=180) for point in range(1, control_points[beams.index(beam)])]
                    gantry.append(gantry_data)
                return gantry
            else: return None
        else:
            mu = []     
            for beam in beams: 
                mu_data = [dicom.BeamSequence[beam].ControlPointSequence[point].CumulativeMetersetWeight for point in range(1, control_points[beams.index(beam)])]
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
                leaf_num: int, y1: float, y2: float
            ) -> bool:
            """Return a boolean specifying if the given leaf is under one of the y jaws.

            Parameters
            ----------
            is_hdmlc: bool
            leaf_num: int
            y1 : float
            y2 : float
            """
            outer_leaf_thickness = 10  # mm
            inner_leaf_thickness = 5
            mlc_position = 0
            if is_hdmlc:
                outer_leaf_thickness /= 2
                inner_leaf_thickness /= 2
                mlc_position = 100
            for leaf in range(1, leaf_num + 1):
                if 10 >= leaf or leaf >= 110:
                    mlc_position += outer_leaf_thickness
                elif 50 >= leaf or leaf >= 70:
                    mlc_position += inner_leaf_thickness
                else:  # between 50 and 70
                    mlc_position += outer_leaf_thickness

            y2_position = y2 * 10 + 200
            y1_position = 200 - y1 * 10
            if 10 >= leaf_num or leaf_num >= 110:
                thickness = outer_leaf_thickness
            elif 50 >= leaf_num or leaf_num >= 70:
                thickness = inner_leaf_thickness
            else:  # between 50 and 70
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

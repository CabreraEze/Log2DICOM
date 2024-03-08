"""Gamma_Bakai & Snapshot to control point mapping functions"""
from pylinac import TrajectoryLog
from scipy.ndimage import sobel

from bisect import bisect_right
import numpy as np
import copy
import argue


def gamma_Bakai(
            reference: np.ndarray,
            comparison: np.ndarray,
            resolution: float = 0.1,
            threshold: float = 0.1,
            distTA: float = 1,
            doseTA: float = 1,
            report: bool = False
    ) -> np.ndarray:

    ref_img = copy.deepcopy(reference)
    ref_img /= ref_img.max()
    comp_img = copy.deepcopy(comparison)
    comp_img /= comp_img.max()

    dpi = 25.4 / resolution
    dpmm = dpi / 25.4

    ref_img[ref_img < threshold * np.max(ref_img)] = np.NaN     # invalidate dose values below threshold so gamma doesn't calculate over it
    distTA_pixels = dpmm * distTA   # convert distance value from mm to pixels

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
        
        
        def interpol_avoidance(indices):
            interpol_indices = indices
            diff = np.diff(indices)

            if 0 in diff: 
                mean_diff = int(np.mean(diff[diff != 0]))

                zero_index = [i for i,j in enumerate(diff) if j==0]
                for z_id in zero_index:
                    interpol_indices[z_id] = interpol_indices[z_id-1] + mean_diff
                
            return interpol_indices


        indices = [0,]
        for point_meterset in meterset:
            diff = dyna_meterset - point_meterset
            indices.append(bisect_right(diff,0)-1)
        return interpol_avoidance(indices)


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
        # else:
        #     raise popAdvice('¨kind¨ doesnt math any aceptable case')
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
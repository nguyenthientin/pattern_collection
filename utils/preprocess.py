"""
 @author Tin T. Nguyen
 @email t.t.nguyen-3@tudelft.nl
 @create date 2020-11-03
 @modify date 2020-11-12
 @desc [description]
""" 

import warnings
import numpy as np

from constants import SPEED_NaN

def preprocess_speed(speed):
    speed = speed.astype(np.float)
    speed[speed == SPEED_NaN] = np.nan
    isnan_speed = np.isnan(speed)
    if np.all(isnan_speed):
        return []
    # if np.any(isnan_speed):
    #     warnings.warn('NaN(s) found in speed. Replace those with maximum speed.'
    #                    + 'The result might be unreliable')
    max_speed = np.max(speed, initial=-1, where=~isnan_speed)
    speed[isnan_speed] = max_speed
    return speed
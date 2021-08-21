import numpy as np
import pandas as pd


def get_vaccination_rate(stock, horizon_time):
    vaccination_rate = -1.0 * np.log(1.0 - stock) / horizon_time
    return vaccination_rate

# -*- coding: utf-8 -*-
# based on python 3.11
# author haodonLiu

import time
import numpy as np
from functools import wraps


def timer(func):
    @wraps(func)
    def wrapper(*arg, **karg):
        start = time.perf_counter()
        be_dec = func(*arg, **karg)
        duration = time.perf_counter() - start
        print(f'[{wrapper.__name__}] took {duration * 1000} ms')
        return be_dec
    return wrapper


def gui1(item):  # zero scale归一算法
    theta = np.std(np.array(item))
    ave = np.average(np.array(item))
    new_item = []
    for i in item:
        i = (i-ave)/theta
        new_item.append(i)
    return np.array(new_item)

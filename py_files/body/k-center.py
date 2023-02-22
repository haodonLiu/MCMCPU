from dic import set
from pickle import load
import numpy as np
import methods

num = set["num"]

dfin = load(open(set['root']+r"\value.stfi", "rb"))
points = list(zip(dfin['rmsd'], dfin['rg2']))


def get_init_center(k):
    dict = {}
    rn = round(np.random.uniform(0, num-1))
    init_cen = points[rn]
    centerlist = [init_cen]
    dis_lis = []

    for i in centerlist:
        for k in points:
            dis_lis = []
            dis = methods.dis2d(i, k)
            dis_lis.append(dis)
            dict[i] = dis_lis
        s_dis_lis = sorted(dis_lis)
        cen = dict[s_dis_lis[-1]]
        centerlist.append(cen)

        if len(centerlist) == k:
            return centerlist

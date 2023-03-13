from pickle import load
import math_methods



def k_center(k,programe):
    dfin = load(open(programe.path, "rb"))
    points = list(zip(dfin['rmsd'], dfin['rg2']))
    dict = {}
    rn = round(np.random.uniform(0, programe.total_num-1))
    init_cen = points[rn]
    centerlist = [init_cen]
    dis_lis = []

    for i in centerlist:
        for k in points:
            dis_lis = []
            dis = math_methods.dis2d(i, k)
            dis_lis.append(dis)
            dict[i] = dis_lis
        s_dis_lis = sorted(dis_lis)
        cen = dict[s_dis_lis[-1]]
        centerlist.append(cen)

        if len(centerlist) == k:
            return centerlist

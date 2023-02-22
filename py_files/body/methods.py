
def rmsd(df_mark, df, mode=all):
    if mode == 'CA':
        for i in [df['atom'], df_mark['atom']]:
            if "CA" in i:
                loclis = list(zip(i))
    total = 0
    for n in range(1, len(df_mark)-1):
        x1, y1, z1 = df_mark["x_axe"][n], df_mark["y_axe"][n], df_mark["z_axe"][n]
        x2, y2, z2 = df["x_axe"][n], df["y_axe"][n], df["z_axe"][n]
        total += (x1-x2)**2+(y1-y2)**2 + (z1-z2)**2
    rmsd_value = round((total/len(df_mark))**0.5, 5)
    return rmsd_value


def get_mass_center(df):
    sum_weight = 0.0
    sum_loc = [0, 0, 0]
    n = m = -1
    mass_center = []
    dic = {'H': 1.0079, 'C': 12.011, 'N': 14.007,
           'O': 15.999, 'P': 30.974, 'S': 30.06}
    for (x, y, z) in zip(df['x_axe'], df['y_axe'], df['z_axe']):
        n += 1
        for i in dic.keys():
            if i in df['atom'][n]:
                sum_loc[0] += x*dic[i]
                sum_loc[1] += y*dic[i]
                sum_loc[2] += z*dic[i]
                sum_weight += dic[i]

    for i in range(3):
        weighting = sum_weight
        mass_center.append(sum_loc[i]/weighting)

    mc = mass_center
    w_dis2 = 0
    for (x, y, z) in zip(df['x_axe'], df['y_axe'], df['z_axe']):
        dis2 = (mc[0]-x)**2+(mc[1]-y)**2+(mc[2]-z)**2
        m += 1
        for i in dic.keys():
            if i in df['atom'][m]:
                w_dis2 = dis2*dic[i]
    rg2 = w_dis2/sum_weight

    return mass_center, rg2


def dis2d(a, b):
    x1, y1 = a[0], a[1]
    x2, y2 = b[0], b[1]
    dis2 = (x1-x2)**2+(y1-y2)**2
    dis = dis2**0.5
    return dis

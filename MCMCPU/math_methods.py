import numpy as np

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




def dis2d(a, b):
    aloc = np.array((a[0],a[1]))
    bloc = np.array((b[0],b[1]))
    gap=aloc-bloc
    dis = np.sqrt(gap.dot(gap))
    return dis

def dis3d(a, b):
    aloc = np.array((a[0],a[1],b[2]))
    bloc = np.array((b[0],b[1],b[2]))
    gap=aloc-bloc
    dis = np.sqrt(gap.dot(gap))
    return dis

def gui1(item):  # zero scale归一算法
    theta = np.std(np.array(item))
    ave = np.average(np.array(item))
    new_item = []
    for i in item:
        i = (i-ave)/theta
        new_item.append(i)
    return np.array(new_item)
from matplotlib import pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import time


def elc(resd):  # 判断极性？
    da = str(resd)
    bg = da.find(' ')+1
    end = da.find('h')-1
    res = da[bg:end]
    if res in ['VAL', 'ALA', 'LEU', 'ILE', 'PRO', 'MET', 'TRP', 'TRP', 'PHE']:
        elc = -1
    else:
        elc = 1
    bg = da.find('resseq=')+7
    end = da.find(' icode=')
    res_seq = da[bg:end]
    return elc, res_seq


def vpr(vec):  # 处理vector类型为列表
    vet = str(vec)
    bg = vet.find(' ')
    ed = vet.find('>')
    vet = vet[bg:ed].replace(' ', '')
    vet = vet.split(',')
    vector = []
    for i in vet:
        i = float(i)
        vector.append(i)
    return vector


def vects(atom, cen_mas):  # 坐标提取,输出为列表，并计算到质心的距离
    vec = atom.get_vector()
    vdis = vec-cen_mas
    vec = vpr(vec)
    ldis = vpr(vdis)
    sq_sum = 0
    for i in ldis:
        sq_sum += float(i)**2
    d = float(sq_sum**0.5)
    return vec, d


def sequ(atom):  # 原子序列读取
    res = str(atom.get_parent())
    bg = res.find('resseq')+7
    ed = res.find('icode')
    atom_seq = float(res[bg:ed])
    return atom_seq


def avefu(atom_seq_list, item_list):  # 算每个氨基酸的平均
    sum = 0
    n = 1
    ave_list = []
    for i in range(len(atom_seq_list)):
        if i == 0:
            n += 1
            sum += item_list[i]
        elif atom_seq_list[i-1] == atom_seq_list[i]:
            sum += item_list[i]
            n += 1
        else:
            ave = sum/n
            ave_list.append(ave)
            sum = 0
            n = 1
    return ave_list


def main(infi):  # 主题，获得vec_lis（原子向量/位置），e_list（氨基酸极性标记），dis_ave_list（分子内原子到质心平均距离）,b_factor_list，seq_list

    p = PDBParser()
    structure = p.get_structure("test", infi)
    cen_mas = structure.center_of_mass()
    residues = structure.get_residues()
    atomlist = structure.get_atoms()
    veclis = []
    elist = []
    dis_list = []
    atom_seq_list = []
    b_factor_list = []
    res_seq_list = []

    for atom in atomlist:  # 原子层面
        pos, dis = vects(atom, cen_mas)
        veclis.append(pos)
        b_factor = float(atom.get_bfactor())
        atom_seq = sequ(atom)
        dis_list.append(dis)
        atom_seq_list.append(atom_seq)
        b_factor_list.append(b_factor)

    # print(residues)

    for i in residues:  # 残基层面
        e, rse_seq = elc(i)
        elist.append(e)
        res_seq_list.append(rse_seq)

    dis_ave_list = avefu(atom_seq_list, dis_list)
    b_factor_ave_list = avefu(atom_seq_list, b_factor_list)

    return veclis, elist, dis_ave_list, b_factor_ave_list, atom_seq_list, res_seq_list


def gui1(item):  # zero scale归一算法
    theta = np.std(np.array(item))
    ave = np.average(np.array(item))
    new_item = []
    for i in item:
        i = (i-ave)/theta
        new_item.append(i)
    return np.array(new_item)


t1 = time.time()
vec_lis, e_list, dis_ave_list, b_factor_list, atom_seq_list, res_seq_list = main(
    R'C:\Users\LHD\Desktop\idec 2023 储备\pdb\output\new_0001.pdb')

del e_list[-1::]  # 数据长度对不上

e_list_p = gui1(e_list)  # 归一化
dis_ave_list_p = gui1(dis_ave_list)
b_factor_list_p = gui1(b_factor_list)
t2 = time.time()


def vis(e, x, y, fit=False, table=True, x_name='x', y_nane='y'):
    fig, axlist = plt.subplots(1, 2, sharex=True, sharey=True)
    xmi = min(x)
    ymi = min(y)
    xma = max(x)
    yma = max(y)

    axlist[0].set_xlim(xmi*1.1, xma*1.1)
    axlist[0].set_ylim(ymi*1.1, yma*1.1)

    axlist[0].set_title('nonpolarity')
    axlist[1].set_title('polarity')

    for i in axlist:
        i.set_xlabel(x_name)
        i.set_ylabel(y_nane)

    polist = []
    ipolist = []

    for i in range(len(e)):
        if e[i] > 0:
            polist.append([x[i], y[i], res_seq_list[i]])
        if e[i] < 0:
            ipolist.append([x[i], y[i], res_seq_list[i]])

    for i in polist:
        axlist[0].scatter(i[0], i[1])
    for i in ipolist:
        axlist[1].scatter(i[0], i[1])

    if fit == True:
        def func(x, a, b, c):
            return a * x**2+b*x + c
        popt, pcov = curve_fit(func, x, y)
        yy2 = [func(i, popt[0], popt[1], popt[2]) for i in x]
        plt.plot(x, yy2, 'r-', label='fit')

    if table == True:
        for i in polist:
            axlist[0].text(i[0], i[1], i[2])
        for i in ipolist:
            axlist[1].text(i[0], i[1], i[2])

    plt.show()


vis(e_list_p, dis_ave_list_p, b_factor_list_p,
    table=False, x_name='dis', y_nane='b_fac')

print(t2-t1)

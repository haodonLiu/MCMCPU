# ['atom_num', 'atom', 'res', 'res_num', 'x_axe', 'y_axe', 'z_axe']elements got from pdb
def main(filein, AAs, mode='mg'):
    tdic = {}
    # initialization

    file = open(filein, 'rt')
    content = file.readlines()
    file.close()
    # open files

    if mode[0] == 'g':
        a = 0
    else:
        a = 1

    if mode[1] == 'g':
        b = 0
    else:
        b = 1
    # mode:m=mass/g=geometry,[0]:per AA,[1]:AAs

    for line in content:
        if "ATOM" in line:
            atom_num = int(line[7:11])
            atom = line[77:78]
            res = line[17:20]
            res_num = int(line[22:26].replace(' ', '0'))
            x_axe = float(line[30:38].replace('.', ''))/1000
            y_axe = float(line[39:47].replace('.', ''))/1000
            z_axe = float(line[48:54].replace('.', ''))/1000
            # 以上都是对pdb文件切片获取数据
            if "CA" in line[12:17]:
                atom = "CA"
                # special mark for "CA"

            tdic[atom_num] = [atom_num, atom, res, res_num, x_axe, y_axe, z_axe]

    weight = {'H': 1.0079, 'C': 12.011, 'CA': 12.011,
              'N': 14.007, 'O': 15.999, 'P': 30.974, 'S': 30.06}

    loc_dic = {}
    print(AAs)
    for i in AAs:
        i = int(i)
        nu = 0
        sum_w = 0
        sum_x, sum_y, sum_z = 0, 0, 0

        for n in tdic.values():

            if i in n:
                if a == 1:
                    nu += 1
                    sum_x += n[4]*weight[n[1]]
                    sum_y += n[5]*weight[n[1]]
                    sum_z += n[6]*weight[n[1]]
                    sum_w += weight[n[1]]
                    print(nu)
                else:
                    nu += 1
                    sum_x += n[4]
                    sum_y += n[5]
                    sum_z += n[6]
                    sum_w += weight[n[1]]
            else:
                continue
        if a == 1:
            nu = sum_w

        c_x = round(sum_x/nu, 3)
        c_y = round(sum_y/nu, 3)
        c_z = round(sum_z/nu, 3)
        loc_dic[i] = [c_x, c_y, c_z, round(sum_w, 5)]

    # calculate AA itself center
    t_w = 0
    for i in AAs:
        t_w += loc_dic[i][3]
    #calculate total weight
    sum_x, sum_y, sum_z = 0, 0, 0
    for i in AAs:
        if b == 1:
            sum_x += loc_dic[i][0]*loc_dic[i][3]
            sum_y += loc_dic[i][1]*loc_dic[i][3]
            sum_z += loc_dic[i][2]*loc_dic[i][3]
        else:
            sum_x += loc_dic[i][0]
            sum_y += loc_dic[i][1]
            sum_z += loc_dic[i][2]

    n = len(AAs)

    if b == 1:
        n = t_w

    c_x = sum_x/n
    c_y = sum_y/n
    c_z = sum_z/n
    return ('{:.3f},{:.3f},{:.3f}'.format(c_x, c_y, c_z), mode)


if __name__ == '__main__':
    filein = r'' + \
        input(r'path(whole path without ",you dont need print\,such as:C:xxx\xxx):')
    AAsi = list(input('AA_number,split with ",",such as:76, 183:').split(','))
    AAs = []
    for i in AAsi:
        AAs.append(int(i))
    lis = []
    num = 0
    for i in ['m', 'g']:
        for n in ['m', 'g']:
            mode = i+n
            lis.append(main(filein, AAs, mode))
    fileout = open(filein+'_out.txt', 'wt')
    for i in lis:
        fileout.write('mode={},AAs={},CC={}\n'.format(i[1], AAs, i[0]))
    fileout.close()
    print("finish!output:{}".format(filein+'_out.txt'))

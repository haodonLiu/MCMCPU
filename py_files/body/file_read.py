# -*- coding: utf-8 -*-
# based on python 3.11
# author haodonLiu

import pandas as pd
import pickle
from dic import set
import multiprocessing as mp

# input pdb and transfer into dataframe,save as "stored_file_xxxx(index).stfi"via pickle


class read:
    def readpdb_to_df(self, filein, dc, check):
        # dc hold the total loc of each atom
        col = ['atom_num', 'atom', 'res', 'res_num', 'x_axe', 'y_axe', 'z_axe']
        content = filein.readlines()
        df_raw = []
        dic = {}
        for line in content:
            if "ATOM" in line:
                atom_num = int(line[7:11])
                atom = line[77:78]
                res = line[17:20]
                res_num = int(line[21:26].replace(' ', '0'))
                x_axe = int(line[30:38].replace('.', ''))/1000
                y_axe = float(line[39:47].replace('.', ''))/1000
                z_axe = float(line[48:54].replace('.', ''))/1000
                if "CA" in line[12:17]:
                    atom = "CA"
                df_raw.append(
                    [atom_num, atom, res, res_num, x_axe, y_axe, z_axe])

                dic[atom_num] = [atom_num, atom, res,
                                 res_num, x_axe, y_axe, z_axe]

                if check == False:
                    dc = dic

                if check == True:
                    dc[atom_num][4] += x_axe
                    dc[atom_num][5] += y_axe
                    dc[atom_num][6] += z_axe

                else:
                    continue

        check = True

        filein.close()
        df = pd.DataFrame(df_raw, columns=col)

        return df, dc, check


def read_main():
    pro_df = ''
    num = set["num"]
    dc = {}
    reader = read()
    check = False
    for i in range(num):
        i += 1
        filein = open(set["root"]+r"\output\new_" +
                      "{:0>4}".format(i)+'.pdb', 'rt')
        outfi = open(set["root"]+r"\stored_files\stored_file" +
                     '_'+str(i)+".stfi", "wb")

        pro_df, dc, check = reader.readpdb_to_df(filein, dc, check)

        pickle.dump(pro_df, outfi)  # chucun
        outfi.close()

        if round(100*i/num, 1) in [5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]:
            print("{:.0f}%".format(100*i/num))

    for i in dc.keys():
        for k in range(4, 7):
            dc[i][k] = round(dc[i][k]/num, 3)

    outfi = open(set["root"]+r"\comfo"+'_ave_'+str(num)+".stfi", "wb")
    df_ave = pd.DataFrame(dc, index=['atom_num', 'atom', 'res', 'res_num',
                          'x_axe', 'y_axe', 'z_axe'], columns=[i for i in range(1, len(dc)+1)])
    pickle.dump(df_ave.T, outfi)

    print('final')


if __name__ == '__main__':
    read_main()

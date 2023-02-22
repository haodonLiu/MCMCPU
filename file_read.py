# -*- coding: utf-8 -*-
# based on python 3.11
# author haodongLiu

import pandas as pd
from set import set


# input pdb and transfer into dataframe,save as "stored_file_xxxx(index).stfi"via pickle


def readpdb2df2pic( filein,outpath,sum_dic):
    # dc hold the total loc of each atom
    col = ['atom_num', 'atom', 'res', 'res_num', 'x_axe', 'y_axe', 'z_axe']
    #elements got from pdb
    dic = {}
    #initialization
    content = filein.readlines()
    for line in content:
        if "ATOM" in line:
            atom_num = int(line[7:11])
            atom = line[77:78]
            res = line[17:20]
            res_num = int(line[21:26].replace(' ', '0'))
            x_axe = int(line[30:38].replace('.', ''))/1000
            y_axe = float(line[39:47].replace('.', ''))/1000
            z_axe = float(line[48:54].replace('.', ''))/1000
            #以上都是对pdb文件切片获取数据
            if "CA" in line[12:17]:
                atom = "CA" 
                #special mark for "CA"

            c=0
            if len(sum_dic)<=atom_num and c==0:
                sum_dic[atom_num] = [atom_num, atom, res,res_num, x_axe, y_axe, z_axe]
                if len(sum_dic)==atom_num :
                    c=1
            if len(sum_dic)!=0 and c==1:
                sum_dic[atom_num][4]+=x_axe
                sum_dic[atom_num][5]+=y_axe
                sum_dic[atom_num][6]+=z_axe

            dic[atom_num] = [atom_num, atom, res,res_num, x_axe, y_axe, z_axe]
    df = pd.DataFrame(data=dic,index=col).T
    print(df)
    pd.to_pickle(df,outpath)#将df储存为pickel文件

    return sum_dic

def read_main(num=set['total_num']):
    sum_dic= {}
    sumn=0
    #初始化
   
    for i in set['filelist']:
        
        sumn += 1
        
        filein = open(set["path"]+r"\input"+'\\' +i, 'rt')
        outfi = open(set["path"]+r"\output"+'\\'+i+'.outfile','wb' )
        sum_dic=readpdb2df2pic(filein,outfi,sum_dic)
        filein.close()
        outfi.close()
        
 
     #对input文件夹内所有文件进行遍历,输出每个原子的坐标加和，平均得到平均结构
    
    for i in sum_dic.keys():
            for n in range(4,7):
                sum_dic[i][n]/=num
                sum_dic[i][n]=round(sum_dic[i][n],3)
    ave_dic=sum_dic

    #对坐标平均
    
    outfi_pic = open(set["path"]+r"\information\comformation"+'_ave_'+str(num)+".outfile", "wb")
    outfi_csv = open(set["path"]+r"\information\comformation"+'_ave_'+str(num)+".csv", "wt")
    df_ave = pd.DataFrame(ave_dic, index=['atom_num', 'atom', 'res', 'res_num',
                          'x_axe', 'y_axe', 'z_axe'], columns=[i for i in range(1, len(ave_dic)+1)])
    df_avet=df_ave.T
    df_avet.to_pickle(outfi_pic)
    df_avet.to_csv(outfi_csv)                                                        
    outfi_pic.close()
    outfi_csv.close()
    #将平均数据以pic文件以及csv储存
    print('final')


if __name__ == '__main__':
    read_main(3)

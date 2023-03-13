# -*- coding: utf-8 -*-
# based on python 3.11
# author haodongLiu

import pandas as pd
import pickle as pic

# input pdb and transfer into dataframe,save as "stored_file_xxxx(index).stfi"via pickle


def readpdb2df2pic( filein,outpath,sum_dic):
    # dc hold the total loc of each atom
    col = ['atom_num', 'atom', 'res', 'res_num', 'x_axe', 'y_axe', 'z_axe']
    AA_dic={'UNL':'!','GLY':'G','ALA':'A','VAL':'V','LEU':'L','ILE':'I','MET':'M','SEC':'U','PYL':'O','PRO':'P','SER':'S','CYS':'C','ASN':'N','GLN':'Q','THR':'T','PHE':'F','TYR':'Y','TRP':'W','ASP':'D','GLU':'E','ARG':'R','LYS':'K','HIS':'H'}
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
            aa=AA_dic[res]
            #以上都是对pdb文件切片获取数据
            if "CA" in line[12:17]:
                atom = "CA" 
                #special mark for "CA"

            c=0
            if len(sum_dic)<=atom_num and c==0:
                sum_dic[atom_num] = [atom_num, atom, res,aa,res_num, x_axe, y_axe, z_axe]
                if len(sum_dic)==atom_num :
                    c=1
            if len(sum_dic)!=0 and c==1:
                sum_dic[atom_num][5]+=x_axe
                sum_dic[atom_num][6]+=y_axe
                sum_dic[atom_num][7]+=z_axe

            dic[atom_num] = [atom_num, atom, res,res_num, x_axe, y_axe, z_axe]
    df = pd.DataFrame(data=dic,index=col).T
    pd.to_pickle(df,outpath)#将df储存为pickel文件

    return sum_dic

def read_main(programe):
    num=programe.total_num
    sum_dic= {}
    sumn=0
    #初始化
   
    for i in programe.filelist:
        
        sumn += 1
        
        filein = open(programe.path+r"\input"+'\\' +i, 'rt')
        outfi = open(programe.path+r"\output"+'\\'+i+'.outfile','wb' )
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
    
    outfi_pic = open(programe.path+r"\information\comformation"+'_ave_'+str(num)+".outfile", "wb")
    outfi_csv = open(programe.path+r"\information\comformation"+'_ave_'+str(num)+".csv", "wt")
    df_ave = pd.DataFrame(ave_dic, index=['atom_num', 'atom', 'res','aa', 'res_num',
                          'x_axe', 'y_axe', 'z_axe'], columns=[i for i in range(1, len(ave_dic)+1)])
    df_avet=df_ave.T
    df_avet.to_pickle(outfi_pic)
    df_avet.to_csv(outfi_csv)                                                        
    outfi_pic.close()
    outfi_csv.close()
    #将平均数据以pic文件以及csv储存
    print('function read_main final!')

def get_mass_center(programe):
    read_dic={}
    for filename in programe.filelist:
        filein = open(programe.path+r"\output"+'\\'+filename+'.outfile','rb')
        df=pic.load(filein)
        filein.close()

        sum_weight = 0.0
        sum_loc = [0, 0, 0]
        n , m = 0,0
        mass_center = []
        
        dic = {'H': 1.0079, 'C': 12.011, 'N': 14.007,
            'O': 15.999, 'P': 30.974, 'S': 30.06,'CA': 12.011}
        
        for (x, y, z) in zip(df['x_axe'], df['y_axe'], df['z_axe']):
            n+=1
            sum_loc[0] += x*dic[df['atom'][n]]
            sum_loc[1] += y*dic[df['atom'][n]]
            sum_loc[2] += z*dic[df['atom'][n]]
            sum_weight += dic[df['atom'][n]]

        for i in range(3):
            mass_center.append(sum_loc[i]/sum_weight)

        mc = mass_center
        w_dis2 = 0
        
        for (x, y, z) in zip(df['x_axe'], df['y_axe'], df['z_axe']):
            dis2 = (mc[0]-x)**2+(mc[1]-y)**2+(mc[2]-z)**2
            m += 1
            for i in dic.keys():
                if i in df['atom'][m]:
                    w_dis2 = dis2*dic[i]
        
        rg2 = w_dis2/sum_weight
        read_dic[filename]=[mass_center,rg2]
        
    
    data=pd.DataFrame(read_dic,index=['mass_center','rg2']).T
    outfile = open(programe.path+r"\information"+'\\'+'mass_center_rg2.csv','wt' )
    data.to_csv(outfile)
    outfile.close()
    print('function get_mass_center final!')



def read_gro2df(path):#update 3/13
    with open(path) as filein:
            content = filein.readlines()
    
    # dc hold the total loc of each atom
    
    col = ['atom','atom_num', 'res', 'res_num', 'x_axe', 'y_axe', 'z_axe']

    #elements got from gro
    
    dic = {}

    for line in content:
        if 'MOL' in line:
            atom_num = int(line[15:20])
            atom = line[10:15].replace(' ', '')
            atom=''.join(filter(str.isalpha, atom))
            res = line[5:10].replace(' ', '')
            res_num = int(line[0:5].replace(' ', '0'))
            x_axe = float(line[20:28].replace('.', ''))/1000
            y_axe = float(line[28:36].replace('.', ''))/1000
            z_axe = float(line[36:44].replace('.', ''))/1000
            #以上都是对gro文件切片获取数据

            dic[atom_num] = [atom,atom_num, res,res_num, x_axe, y_axe, z_axe]
    
    df = pd.DataFrame(data=dic,index=col).T
    
    return df


def read_gro_process(programe):
    
    for i in programe.filelist:
        
        filein = open(programe.path+r"\input"+'\\' +i, 'rt')
        outfi = open(programe.path+r"\output"+'\\'+i+'_re.csv','wt' )
        df=read_gro2df(filein)
        
        sum_weight = 0.0
        sum_loc = [0, 0, 0]
        n = 0
        mass_center = []
        dic = {'H': 1.0079, 'C': 12.011, 'N': 14.007,
            'O': 15.999, 'P': 30.974, 'S': 30.06,'CA': 12.011}
        
        for (x, y, z) in zip(df['x_axe'], df['y_axe'], df['z_axe']):
            n+=1
            sum_loc[0] += x*dic[df['atom'][n]]
            sum_loc[1] += y*dic[df['atom'][n]]
            sum_loc[2] += z*dic[df['atom'][n]]
            sum_weight += dic[df['atom'][n]]

        for i in range(3):
            mass_center.append(sum_loc[i]/sum_weight)
        
        for (x, y, z) in zip(df['x_axe'], df['y_axe'], df['z_axe']):
            df['x_axe'][n]-=mass_center[0]
            df['y_axe'][n]-=mass_center[1]
            df['z_axe'][n]-=mass_center[2]
            n-=1
        df.to_csv(outfi)
        filein.close()
        outfi.close()
    
def top_read(content:list):
    n=-1
    atom_dic={}
    atom=0
    bond=0
    atomname=''
    atom_num=0
    diend=0
    bend=0
    diherdral=0
    for line in content:
        n+=1
        if '[ moleculetype ] ' in line:
            atomname=''.join(filter(str.isalpha,content[n+2]))
            atomnum=''.join(filter(str.isnumeric,content[n+2]))
        if '[ atoms ]' in line:
            atom=n
        if '[ bonds ]' in line:
            bond=n
        if bond>0 and bend==0 and len(line)==1:
            bend=0
            
        if ';   ai    aj    ak    al    funct definition' in line:
            diherdral=n
        if '[ position_restraints ]' in line:
            diend=n

        
        
    for line in content[atom+1:bond-1]:
        atom_num = int(line[3:6])
        atom = line[6]
        res = line[12:15]
        res_num = ''.join(filter(str.isnumeric,line[7:12]))
        atom_dic[atom_num]=[atom,atom_num,res,res_num]
    
    link_dic={}
    
    for i in range(atom_num):
        link_dic[i]=[]
    
    for line in content[bond+1:bend]:
        lis=line.split(' ')        
        lis=list(set(lis))
        lis.remove('')
        ac1=atom_dic[lis[0]][0]+str(atom_dic[lis[0]][0])
        ac2=atom_dic[lis[1]][0]+str(atom_dic[lis[1]][0])
        ac2.lower()
        link_dic[ac1].append(ac2)
        ac2=atom_dic[lis[0]][0]+str(atom_dic[lis[0]][0])
        ac1=atom_dic[lis[1]][0]+str(atom_dic[lis[1]][0])
        ac2.lower()
        link_dic[ac1].append(ac2)
    
    dih_lis=[]
    for line in content[diherdral+1:diend-1]:
        if [ 'dihedrals' ] in line:
            continue
        lis=line.split(' ')        
        lis=list(set(lis))
        lis.remove('')
        ac1=atom_dic[lis[0]][0]+str(atom_dic[lis[0]][0])
        ac2=atom_dic[lis[1]][0]+str(atom_dic[lis[1]][0])

    return link_dic
            
        
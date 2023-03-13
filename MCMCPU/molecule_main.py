import os
import read_gro
import move


class programe():
    def __init__(self):
        self.path= r"C:\Users\LHD\Desktop\IDEC 2023\Conformation\molecule"
        self.times=10
        self.delta=[x/10 for x in range(21)]#[0.0,0.1__2.0]
        self.concat=True
        self.creat_csv=False
        self.link_h_ignore=False
        #一些基本参数path文件路径,'0':1创建字符：数字对应关系
        
        for i in ['\\input','\\output','\\information']:
            check=os.path.exists(self.path+i)
            if check!=True:
                os.mkdir(self.path+i)
        #检验当前目录下是否有input/output/information文件夹

        filelist=os.listdir(self.path+'\\input')
        #遍历所有文件
        for i in filelist:
            if '.gro' in i :
                continue
            else :
             filelist.remove(i)
        self.standard=filelist[0]

        #只读取pdb文件

        total_num=len(filelist)
        self.total_num=total_num
        self.filelist=filelist
        #整理后传参


if __name__=='__main__':
    pro=programe()

    move.diheral(pro,False)



    

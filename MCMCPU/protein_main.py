import os
import protein_pdb_read


class programe():
    def __init__(self) -> None:
        self.path= r"C:\Users\LHD\Desktop\IDEC 2023\Conformation\protein\30_test"
        self.times=1
        self.name="30_test"
        #一些基本参数path文件路径，times聚类次数，nanmes对项目命名
        for i in ['\\input','\\output','\\information']:
            check=os.path.exists(self.path+i)
            if check!=True:
                os.mkdir(self.path+i)
        #检验当前目录下是否有input/output/information文件夹

        filelist=os.listdir(self.path+'\\input')
        #遍历所有文件
        for i in filelist:
            if '.pdb' in i :
                continue
            else :
             filelist.remove(i)
        #只读取pdb文件

        total_num=len(filelist)
        self.total_num=total_num
        self.filelist=filelist
        #整理后传参

        for i in [self.path,self.name,self.total_num,self.times,self.filelist]:
            print(i)


if __name__=='__main__':
    pro=programe()
    protein_pdb_read.read_main(pro)
    protein_pdb_read.get_mass_center(pro)

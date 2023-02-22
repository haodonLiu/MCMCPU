import os

set = {
    'path': r"C:\Users\LHD\Desktop\IDEC 2023\Conformation\30_test",
    'times':1,
    'name':"30_test"
}
#一些基本参数path文件路径，times聚类次数，nanmes对项目命名

for i in ['\\input','\\output','\\information']:
    check=os.path.exists(set["path"]+i)
    if check!=True:
        os.mkdir(set["path"]+i)
#检验当前目录下是否有input/output/information文件夹

filelist=os.listdir(set["path"]+'\\input')
#遍历所有文件
for i in filelist:
    if '.pdb' in i :
        continue
    else :
        filelist.remove(i)
#只保留pdb文件

total_num=len(filelist)
set['total_num']=total_num
set["filelist"]=filelist
#整理后传参

for i in ['path','times','name','total_num']:
    print(i,':',set[i])

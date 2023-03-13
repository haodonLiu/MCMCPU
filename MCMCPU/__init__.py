import os
import pandas as pd
import numpy as np
import mcmcpu_read as read
import link

class project:
    def creat_project(self):
        project_path=input()
    def __init__(self,project_path) -> None:
        self.path=os.path.join(r'',project_path)
        for i in ['input','output','information']:
            check=os.path.exists(os.path.join(self.path,i))
            if check!=True:
                os.makedirs(os.path.join(self.path,i))




class molecule(project):
    def __init__(self,molename,filetype,project_path) -> None:
        project.__init__(self,project_path)
        self.mol_name=str(molename)
        self.type=filetype
        self.mol_path=os.path.join(self.path,'input',self.mol_name+'.'+self.type)
        self.df=read.read_gro2df(self.mol_path)  
                
    def dataframe(self):
        return print(self.df)
    
    def exp_link(self):
        link_dic_exp=link.link_check(self.df,self.mol_name)
        return link_dic_exp



project_path=r'C:\Users\LHD\Desktop\IDEC 2023\Conformation\molecule'
pro=project(project_path)   
mol=molecule('dopa','gro',project_path)
print(mol.exp_link())


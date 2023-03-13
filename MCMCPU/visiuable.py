import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import pickle as pic


def visiuable(programe):
    filein = open(programe.set["path"]+r"\information"+'\\'+'mass_center_rg2.csv','rt')
    dataset=
    filein.close()
   
    dic_file.close()
    rmsd_lis=dataset['rmsd']
    rg2_lis=dataset['rg2']
    plt.scatter(rmsd_lis,rg2_lis)
    plt.show()
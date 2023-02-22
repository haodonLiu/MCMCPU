import seaborn as sns
import matplotlib.pyplot as plt
import pickle as pic
from dic import set


dic_file=open(set['root']+r"\value.stfi","rb")
dataset=pic.load(dic_file)
dic_file.close()
rmsd_lis=dataset['rmsd']
rg2_lis=dataset['rg2']
plt.scatter(rmsd_lis,rg2_lis)
plt.show()
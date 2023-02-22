import methods
from pickle import load
from dic import set
k2 = [([4.21072, 0.04630913], 1955), ([4.90058, 0.04617267], 2361)]
k3 = [([4.02975, 0.04619187], 2208),
      ([4.50905, 0.04617661], 2996), ([5.11949, 0.04639001], 969)]

num = set["num"]


d = {}
for i in [1995, 2361]:
    for k in [2208, 2996, 969]:
        df1 = load(
            open(set['root']+r"\stored_files\stored_file_"+str(i)+".stfi", "rb"))
        df2 = load(
            open(set['root']+r"\stored_files\stored_file_"+str(k)+".stfi", "rb"))
        d[i, k] = methods.rmsd(df1, df2)
print(d)

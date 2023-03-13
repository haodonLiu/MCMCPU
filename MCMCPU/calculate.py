import pickle as pic
import math_methods
import protein_pdb_read as ppr
from set import set
import pandas as pd


def main(num):
    file = open(set['path']+r"\information\comformation_ave_30.outfile", "rb")
    df_mark = pic.load(file)
    file.close()
    #读取平均结构
    rmsd_lis = []
    rg2_lis = []
    mass_center_lis = []

    for i in set['filelist']:

        file = open(set['path']+"\\output\\" +i+".outfile", "rb")
        df = pic.load(file)
        file.close()
        #读取结构
        rmsd = math_methods.rmsd(df_mark, df)
        print(df_mark)
        rmsd_lis.append(rmsd)
        #计算rsmd
        mass_center, rg2 = ppr.get_mass_center(df)
        mass_center_lis.append(mass_center)
        rg2_lis.append(rg2)


    dic_file = open(set['path']+r"\value.stfi", "wb")
    dataset = [rmsd_lis, rg2_lis, mass_center_lis]
    name_lis = ['rmsd', 'rg2', 'mass_of_center']
    col = [i for i in range(1, num+1)]
    dfs = pd.DataFrame(dataset, index=name_lis, columns=col).T
    pic.dump(dfs, dic_file)
    dic_file.close()
    dfs.to_csv(set['path']+r"\information\value.csv")

    print('calculate final!')


if __name__ == '__main__':
    main(set['total_num'])

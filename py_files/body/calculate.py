import pickle as pic
import methods
from dic import set
import pandas as pd


def main(num):
    file = open(set['root']+r"\comfo_ave_3000.stfi", "rb")
    df_mark = pic.load(file)
    file.close()

    rmsd_lis = []
    rg2_lis = []
    mass_center_lis = []

    for i in range(num):

        file = open(set['root']+r"\stored_files\stored_file_" +
                    str(i+1)+".stfi", "rb")
        df = pic.load(file)
        file.close()

        rmsd = methods.rmsd(df_mark, df)
        print(df_mark)
        rmsd_lis.append(rmsd)

        mass_center, rg2 = methods.get_mass_center(df)
        mass_center_lis.append(mass_center)
        rg2_lis.append(rg2)

        if round(100*i/num, 1) in [5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]:
            print("{:.0f}%".format(100*i/num))

    dic_file = open(set['root']+r"\value.stfi", "wb")
    dataset = [rmsd_lis, rg2_lis, mass_center_lis]
    name_lis = ['rmsd', 'rg2', 'mass_of_center']
    col = [i for i in range(1, num+1)]
    dfs = pd.DataFrame(dataset, index=name_lis, columns=col).T
    pic.dump(dfs, dic_file)
    dic_file.close()
    dfs.to_csv(set['root']+r"\value.csv")

    print('final')


if __name__ == '__main__':
    main(set['num'])

import os

import link
import numpy as np
import pandas as pd


def re_read_gro2df(df,times,programe):

    context=[str(times)]
    n=-1

    for i in df['atom1']:
        n+=1
        context.append('{:5}{:<4}{:5}{:5}{:8.3f}{:8.3f}{:8.3f}'.format(df['res_num1'][n],df['res1'][n],str(df['atom1'][n])+str(df['atom_num1'][n]),df['atom_num1'][n],df['x_axe1'][n],df['y_axe1'][n],df['z_axe1'][n]))
    
    if programe.concat:
        n=-1
        for i in df['atom']:
            n+=1
            context.append('{:5}{:<4}{:5}{:5}{:8.3f}{:8.3f}{:8.3f}'.format(df['res_num'][n],df['res'][n],str(df['atom'][n])+str(df['atom_num'][n]),df['atom_num'][n],df['x_axe'][n],df['y_axe'][n],df['z_axe'][n]))
    
    return context

def standard(programe):
    df_standard=[]
    if programe.concat:
        filein = open(programe.path +
                            r"\output"+'\\'+programe.standard+'_re.csv', 'rt')
        df_standard= pd.read_csv(filein)
        del df_standard['Unnamed: 0']
        filein.close()
    return df_standard

def check(path):
    check=os.path.exists(path)
    if check!=True:
        os.makedirs(path)
            

def translational(programe, delta):

    dic={}    
    row = np.random.randint(-10, 11)/10

    df_standard=standard(programe)
    
    for filename in programe.filelist:
        for t in range(programe.times):
            n = -1

            filein = open(programe.path +r"\output"+'\\'+filename+'_re.csv', 'rt')
            df = pd.read_csv(filein)
            filein.close()

            for (x, y, z) in zip(df['x_axe'], df['y_axe'], df['z_axe']):
                n += 1
                dic[n]=[df['atom'][n], df['res'][n], df['res_num'][n],df['x_axe'][n] + row*delta,df['y_axe'][n] + row*delta,df['z_axe'][n]+row*delta,df['atom_num'][n]]
            df=pd.DataFrame(dic,index= ['atom1', 'res1', 'res_num1', 'x_axe1', 'y_axe1', 'z_axe1','atom_num1']).T
            

            outfi_path=programe.path+r"\output"+'\\'+filename+'_re\\move\\translational\\'+'d='+str(delta)+'concat='+str(programe.concat)
            
            check(outfi_path)

            if programe.concat:
                df=pd.concat([df,df_standard],axis=1)
            
            if programe.creat_csv:
                csv_outfi = open(outfi_path+'\\'+str(t)+'.csv', 'wt')
                df.to_csv(csv_outfi)
                csv_outfi.close()

            context=re_read_gro2df(df,t,programe)
            
            outfi=open(outfi_path+'\\'+str(t)+'.gro', 'wt')
            for i in context:
                outfi.write(i+'\n')
            outfi.close()
    
    print('function translational finish!')


def rotate(programe,theta):
    
    df_standard=standard(programe)   
    
    for filename in programe.filelist:
        for t in range(programe.times):
            n = -1
            a, b, c, d ,s1,s2= 0, 0, 0, 0,0,0
            dic={}


            filein = open(programe.path +
                          r"\output"+'\\'+filename+'_re.csv', 'rt')
            df = pd.read_csv(filein)
            filein.close()

            while 1:
                a= np.random.randint(-10, 11)/10
                b= np.random.randint(-10, 11)/10
                c= np.random.randint(-10, 11)/10
                d= np.random.randint(-10, 11)/10
                s1 = round(a**2+b**2,5)
                s2 = round(c**2+d**2,5)
                if  (s1 < 1)*(s2 < 1):
                    break
            

            q0, q1, q2, q3 = a, b, c*((1-s1)/s2)**0.5, d*((1-s1)/s2)**0.5
            A_matrix = [[q0**2+q1**2-q2**2-q3**2, 2*(q1*q2+q0*q3), 2*(q1*q3+q0*q2)], [2*(q1*q2-q0*q3), 
                q0**2-q1**2+q2**2-q3**2, 2*(q2*q3+q0*q1)], [2*(q1*q3+q0*q2), 2*(q2*q3-q0*q1), q0**2-q1**2-q2**2+q3**2]]
            A_matrix=np.array(A_matrix)
            #creat martix


            for (x, y, z) in zip(df['x_axe'], df['y_axe'], df['z_axe']):
                n += 1
                location = np.array([[df['x_axe'][n]], [df['y_axe'][n]], [df['z_axe'][n]]])
                new_location=location*A_matrix
                dic[n]=[df['atom'][n], df['res'][n], df['res_num'][n],new_location[0,0],new_location[1,0],new_location[2,0],df['atom_num'][n]]
            
            df=pd.DataFrame(dic,index= ['atom1', 'res1', 'res_num1', 'x_axe1', 'y_axe1', 'z_axe1','atom_num1']).T
            
            outfi_path=programe.path+r"\output"+'\\'+filename +'_re\\move\\rotation\\'+'theta='+str(theta)+'comcat='+str(programe.concat)
            
            check(outfi_path)
            
            if programe.concat:
                df=pd.concat([df,df_standard],axis=1)
            
            if programe.creat_csv:
                csv_outfi = open(outfi_path+'\\'+str(t)+'.csv', 'wt')
                df.to_csv(csv_outfi)
                csv_outfi.close()

            context=re_read_gro2df(df,t,programe)
            
            outfi=open(outfi_path+'\\'+str(t)+'.gro', 'wt')
            for i in context:
                outfi.write(i+'\n')
            outfi.close()
    
    print('function rotate finish!')

def diheral(programe,h_ignore:bool):
    
    for filename in programe.filelist:
        link_dic=link.link_check(programe,filename)#use bond length to assume link

        chose_dic=link_dic.copy()#{'O10': ['h11', 'c5'], 'O12': ['h13', 'c6'], 'C14': ['h16', 'h15', 'c17', 'c3'], 'C17': ['h18', 'h19', 'n20', 'c14'], 'N20': ['h22', 'h21']}

        for i in link_dic:#remove invariable bond
            if 'ph' in link_dic[i]:
                del chose_dic[i]

            if 'double' in link_dic[i] :
                del chose_dic[i]
            
            if 'H' in i:
                del chose_dic[i]
        
        variable_diheral_dict={}#{'O10': ['c5'], 'O12': ['c6'], 'C14': ['c3'], 'C17': ['c14']}
        for i in chose_dic:#creat a new list to hold bond with out H atom,A-H isnt a good rorate axe
            
            for n in chose_dic[i]:
                
                if 'h' not in n:
                    variable_diheral_dict[i]=''

                    variable_diheral_dict[i]+=n.upper()
        print(variable_diheral_dict)

        filein = open(programe.path +
                          r"\output"+'\\'+filename+'_re.csv', 'rt')
        df = pd.read_csv(filein)
        filein.close()
        #load the atom location as a dataframe
        
        n=0

        new_lis_n=[]#['10', '5', '12', '6', '14', '3', '17', '14']
        new_lis_an=[]#['O10', 'C5', 'O12', 'C6', 'C14', 'C3', 'C17', 'C14']
        tier_1=link_dic.copy()

        for i in variable_diheral_dict.items():#copy into two list for futher import
            for k in i: 
                new_lis_an.append(k)
                k="".join(filter(str.isdigit, k))
                new_lis_n.append(k)

        for i in range(0,len(new_lis_n),2):
            
            k=int(new_lis_n[i])-1
            an2=new_lis_an[i]
            an3=new_lis_an[i+1]

            if h_ignore:#ususally the atom H shows a low contribution to the system energy
                for i in [an2,an3]:
                    lin_lis=[]
                    for n in tier_1[i]:
                        if 'h' not in n:
                            lin_lis.append(n)
                    tier_1[i]=lin_lis

            if len(tier_1[an2])==1:
                an1=0
                an4=tier_1[an3][1].upper()

            if len(tier_1[an3])==1:
                an1=tier_1[an2][1].upper()
                an4=0

            else:
                an1=tier_1[an2][1].upper()
                an4=tier_1[an3][1].upper()
            
            tier_2=[]


            a=np.array([df['x_axe'][k-1], df['y_axe'][k-1], df['z_axe'][k-1]])
            b=np.array([df['x_axe'][k], df['y_axe'][k], df['z_axe'][k]])
            ab=b-a

            
            
            

            print(an2,tier_1[an2],an3,tier_1[an3])

        

                          


            






        np.random.randint



import pandas as pd
import math_methods as mathm


def link_check(df,filename):

    index=[]
    for n in range(len(df)):
        n=n+1
        index.append(str(df['atom'][n]).lower()+str(df['atom_num'][n]))
    index.reverse()
    df_dis=pd.DataFrame(index=index)

    n=0
    lis=list(zip(df['x_axe'], df['y_axe'], df['z_axe']))

    for (x, y, z) in zip(df['x_axe'], df['y_axe'], df['z_axe']):
        a=lis[n]
        n+=1
        x=0
        dis_list=[]
        for i in df['x_axe']:
            b=lis[x]
            num=mathm.dis3d(a,b)
            dis_list.append(round(num,5))
            x+=1
        dis_list.reverse()
        df_dis[str(df['atom'][n])+str(df['atom_num'][n])]=dis_list
        
    

    link_dic={}
    
    h=0
    c=0
    n=0
    for i in df_dis.columns:
        if 'H' in i:
            h+=1
        if 'C' in i :
            c+=1
        if 'N' in i :
            n+=1
    pi=(2*c+2-h-n)/2
    
    
    for i in df_dis.columns:
        s=df_dis[i].sort_values().index[1:6]
        sv=df_dis[i].sort_values()   
        atom=[]
        
        for n in s:
            atom.append(n)
        
        if 'H' in i:
            
            while 'h' in atom[0]:   
                atom.remove(atom[0])
            
            link_dic[i]=atom[0]

        if 'O' in i :
            if pi>0 and sv[1]>0.2 and 'c' in s[0]: 
                link_dic[i]=atom[0].append('double')

            else:
                link_dic[i]=atom[0:2] 
        
        if 'N' in i :
            link_dic[i]=atom[0:2]

        if 'C' in i :
            
            if 'H' in s[0] and sv[0]>1.2:
                atom.remove(atom[0])

            elif pi>0 and 'c' in s[1]:

                if 'c' in s[2] and sv[2]<0.16:
                    st=atom[0:3]
                    st.append('ben')
                    
                    link_dic[i]=st
                    

                else:
                    st=atom[0:3]
                    st.append('double')
                    link_dic[i]=st

            
            else:
                link_dic[i]=atom[0:4]


    return link_dic

def link_process(link_dic):
    chose_dic=link_dic.copy()#{'O10': ['h11', 'c5'], 'O12': ['h13', 'c6'], 'C14': ['h16', 'h15', 'c17', 'c3'], 'C17': ['h18', 'h19', 'n20', 'c14'], 'N20': ['h22', 'h21']}
    for i in link_dic:#remove invariable bond
        
        if 'ben' in link_dic[i]:
           del chose_dic[i]

        if 'double' in link_dic[i] :
            del chose_dic[i]
            
        if 'H' in i:
            del chose_dic[i]
    
    pro_link_dic=chose_dic
    
    var_diheral_dic={}

    for i in chose_dic:#creat a new list to hold bond with out H atom,A-H isnt a good rorate axe
                
        for n in chose_dic[i]:
            
            if 'h' not in n:
                var_diheral_dic[i]=''
                var_diheral_dic[i]+=n.upper()
    
    
    return pro_link_dic,var_diheral_dic
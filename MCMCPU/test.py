import numpy as np
a=np.array([1,2,5],dtype=np.float64)
b=np.array([2,3,4],dtype=np.float64)
c=np.array([3,4,5],dtype=np.float64)
d=np.array([9,5,5],dtype=np.float64)
def diheral_calculate(a,b,c,d):
    ab=a-b
    bc=b-c
    cd=c-d
    i1=np.cross(ab,bc)
    i2=np.cross(bc,cd)
    _i1_=np.sqrt(i1.dot(i1))
    _i2_=np.sqrt(i2.dot(i2))
    cos=i1.dot(i2)/(_i1_*_i2_)
    angle=(np.pi-np.arccos(cos))*360/(2*np.pi)
    return round(angle,3)

def re_diheral_calculate(diangle,a,b,c,d,delta):
    cos=-np.cos(diangle*(2*np.pi)/360)
    ab=a-b
    bc=b-c
    cd=c-d
    i1=np.cross(ab,bc)
    _i1_=np.sqrt(i1.dot(i1))
    m=round(cos*_i1_,5)
   
    while 1:
        d=(delta+1)*d
        i2=np.cross(bc,cd)
        _i2_=np.sqrt(i2.dot(i2))
        func=round(i1.dot(i2)/_i2_)
        if func==m:
            print(d)
            break

    





delta=0.01
diangle=diheral_calculate(a,b,c,d)
re_diheral_calculate(diangle+1,a,b,c,d,delta)

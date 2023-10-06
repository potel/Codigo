import numpy as np
import matplotlib.pyplot as plt
import sys
import cmath
from matplotlib import rc,rcParams
rc('text',usetex=True)
rc('font',**{'family':'serif','serif':['Computer Modern'],'size':13})
s1="dsdE.txt"
with open(s1) as f1:
    lines = f1.readlines()
    x1 = [float(line.split()[1]) for line in lines]
    y1 = [float(line.split()[2]) for line in lines]
    y2 = [float(line.split()[4]) for line in lines]
    y3 = [float(line.split()[6]) for line in lines]
    y4 = [float(line.split()[8]) for line in lines]
    y5 = [float(line.split()[10]) for line in lines]
    y6 = [float(line.split()[12]) for line in lines]
    y7 = [float(line.split()[14]) for line in lines]
    y8 = [float(line.split()[16]) for line in lines]
    y9 = [float(line.split()[18]) for line in lines]
    y10 = [float(line.split()[20]) for line in lines]
    y11 = [float(line.split()[22]) for line in lines]
    y12 = [float(line.split()[21]) for line in lines]
    y13 = [float(line.split()[23]) for line in lines]
    y14 = [float(line.split()[24]) for line in lines]
#with open(s2) as f1:
#    lines = f1.readlines()
#    x2 = [line.split()[0] for line in lines]
#    y2 = [line.split()[1] for line in lines]
#with open(s3) as f1:
#    lines = f1.readlines()
#    x3 = [line.split()[0] for line in lines]
#    y3 = [line.split()[1] for line in lines]
#with open(s4) as f1:
#    lines = f1.readlines()
#    x4 = [line.split()[0] for line in lines]
#    y4 = [line.split()[1] for line in lines]

#yy=[0. for i in x1]
#count=-1
#for n in x1:
#    yy[count]=4*float(y1[count])
#    count=count+1    
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.set_title("$^{84}$Se $(d,p)$, $E_d=90$ MeV",fontsize=23)
ax1.set_xlabel('$E_x$ (MeV)',fontsize=23)
ax1.set_ylabel('$d\sigma/dE$ (mb/MeV)',fontsize=23)

ax1.plot(x1,y1, c='k',linestyle='--',lw=2,label="NEB $\ell=0$")
ax1.plot(x1,y2, c='b',lw=2,label="NEB $\ell=1$")
ax1.plot(x1,y3, c='g',lw=2,label="NEB $\ell=2$")
ax1.plot(x1,y4, c='r',lw=2,label="NEB $\ell=3$")
ax1.plot(x1,y5, c='k',lw=2,label="NEB $\ell=4$")
ax1.plot(x1,y6, c='r',linestyle='--',lw=2,label="NEB $\ell=5$")
ax1.plot(x1,y7, c='b',linestyle='--',lw=2,label="NEB $\ell=6$")
ax1.plot(x1,y8, c='g',linestyle='--',lw=2,label="NEB $\ell=7$")
ax1.plot(x1,y9, c='r',linestyle='-.',lw=2,label="NEB $\ell=8$")
#ax1.plot(x1,y14, c='k',linestyle='-',lw=3,label="total")
#ax1.plot(x1,y11, c='b',linestyle='-',lw=3,label="total NEB")
#ax1.plot(x1,y13, c='r',linestyle='-',lw=3,label="EB")
#ax1.plot(x1,y10, c='r',lw=2,label="NEB $\ell=9$")
#ax1.plot(x1,y11, c='k',lw=4,label="total NEB")
major_ticksx = np.arange(-10, 10,0.4)
minor_ticksx = np.arange(-10, 90,0.2)
major_ticksy = np.arange(0, 2, 0.2)
minor_ticksy = np.arange(0, 2, 0.1)
ax1.set_xticks(major_ticksx)
ax1.set_xticks(minor_ticksx, minor=True)
ax1.set_yticks(major_ticksy)
ax1.set_yticks(minor_ticksy, minor=True)
ax1.tick_params(axis='both', which='major', labelsize=20,width=1,length=6)
ax1.tick_params(axis='both', which='minor',width=1,length=3)



ax1.set_xlim(-2,2)
ax1.set_ylim(0,1.2)
leg = ax1.legend()
plt.tight_layout()
fig.savefig("dsdE_L.pdf", bbox_inches='tight')
plt.show()

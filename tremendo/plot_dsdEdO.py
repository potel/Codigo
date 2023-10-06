import numpy as np
import matplotlib.pyplot as plt
import sys
import cmath
from matplotlib import rc,rcParams
rc('text',usetex=True)
rc('font',**{'family':'serif','serif':['Computer Modern'],'size':13})
s1="dsdEdO.txt"
s2="dsdEdO_high_9MeV.txt"
s3="dsdEdO_high_8MeV.txt"
#s2="cross_15.txt"
#s3="cross_30.txt"
#s4="../data.csv"
with open(s1) as f1:
    lines = f1.readlines()
    x1 = [float(line.split()[0]) for line in lines]
    y1 = [float(line.split()[1]) for line in lines]
with open(s2) as f1:
    lines = f1.readlines()
    x2 = [float(line.split()[0]) for line in lines]
    y2 = [float(line.split()[1]) for line in lines]
with open(s3) as f1:
    lines = f1.readlines()
    x3 = [float(line.split()[0]) for line in lines]
    y3 = [float(line.split()[1]) for line in lines]
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
ax1.set_yscale('log')
ax1.set_title("$^{157}$Gd $(d,p)$, $E_x=-0.5$ MeV",fontsize=23)
ax1.set_xlabel('$\\theta_{CM}$',fontsize=23)
ax1.set_ylabel('$d\sigma/d\Omega dE$ (mb/MeV sr)',fontsize=23)

ax1.plot(x1,y1, c='k',lw=4,label="$E_d=9$ MeV")
#ax1.plot(x2,y2, c='r',lw=4,label="$E_d=9$ MeV")
#ax1.plot(x3,y3, c='b',lw=4,label="$E_d=8$ MeV")
ax1.tick_params(axis='both', which='major', labelsize=20,width=1,length=6)
ax1.tick_params(axis='both', which='minor',width=1,length=3)
major_ticksx = np.arange(0, 181, 20)
minor_ticksx = np.arange(0, 181, 10)
#major_ticksy = np.arange(0, 100, 1)
#minor_ticksy = np.arange(0, 251, 0.2)
ax1.set_xticks(major_ticksx)
ax1.set_xticks(minor_ticksx, minor=True)

ax1.set_xlim(0,180)
#ax1.set_ylim(0,80)
#ax1.plot(x3,y3, c='b',lw=4,label="DWBA, $R_{max}$=30 fm")
#ax1.plot(x4,y4, c='k',marker='o',lw=0,markersize=12,label="data")
leg = ax1.legend()
plt.tight_layout()
fig.savefig("dsdEdO-9MeV-v2.pdf", bbox_inches='tight')
plt.show()

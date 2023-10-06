#!/usr/bin/python3
import numpy as np
import math
import matplotlib.pyplot as plt
import sys
from matplotlib import rc,rcParams
rc('text',usetex=True)
rc('font',**{'family':'serif','serif':['Computer Modern'],'size':13})
def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False
count=0
with open(sys.argv[1]) as f1:
    lines = f1.readlines()
x0=[]
x1=[]
x2=[]
x3=[]
x4=[]
x5=[]
x6=[]
x7=[]
for i in range(len(lines)):
    if (is_number(lines[i].split()[0])):
        x0.append(float(lines[i].split()[int(sys.argv[2])]))
        x1.append(float(lines[i].split()[int(sys.argv[3])]))
#        x2.append(float(lines[i].split()[int(sys.argv[4])]))
#        x3.append(float(lines[i].split()[int(sys.argv[5])]))
    else:
        print(lines[i])
with open(sys.argv[4]) as f2:
    lines = f2.readlines()
for i in range(len(lines)):
    if (is_number(lines[i].split()[0])):
        x2.append(float(lines[i].split()[int(sys.argv[5])]))
        x3.append(float(lines[i].split()[int(sys.argv[6])]))
    else:
        print(lines[i])
with open(sys.argv[7]) as f3:
    lines = f3.readlines()
for i in range(len(lines)):
    if (is_number(lines[i].split()[0])):
        x4.append(float(lines[i].split()[int(sys.argv[8])]))
        x5.append(float(lines[i].split()[int(sys.argv[9])]))
    else:
        print(lines[i])
with open(sys.argv[10]) as f4:
    lines = f4.readlines()
for i in range(len(lines)):
    if (is_number(lines[i].split()[0])):
        x6.append(float(lines[i].split()[int(sys.argv[11])]))
        x7.append(float(lines[i].split()[int(sys.argv[12])]))
    else:
        print(lines[i])

fig = plt.figure()
ax1 = fig.add_subplot(111)
#ax1.set_xscale('log')
#ax1.set_yscale('log')
xx3=[x3[i]*-6e3 for i in range(len(x3))]
xx1=[x1[i]*-1.5e4 for i in range(len(x1))]
xx5=[x5[i]*-5e4 for i in range(len(x5))]
xx7=[x7[i]*6e5 for i in range(len(x7))]
ax1.set_title('$^{60}$Ni($^{116}$Sn,$^{114}$Sn)$^{62}$Ni, $E_{cm}=154.26$ MeV',fontsize=20)
ax1.set_xlabel('$r_{Aa}$ (fm)',fontsize=20)
ax1.set_ylabel('arbitrary units',fontsize=20)
ax1.plot(x0,xx1,lw=2,c='k',label='$Q_{2n}=1.307$ MeV')
ax1.plot(x2,xx3,lw=2,c='b',label='$Q_{2n}=0.1$ MeV')
ax1.plot(x4,xx5,lw=2,c='g',label='$Q_{2n}=-8$ MeV')
ax1.plot(x6,xx7,lw=2,c='r',label='$Q_{2n}=-20$ MeV')
#ax1.plot(x2,xx3,lw=3,c='r',label='gauge phase effect ON')
#ax1.plot(x0,x1,lw=3,c='k',label='$\gamma$ probability')
#ax1.plot(x0,x2,lw=3,c='b',label='$\\nu$ probability')
#ax1.plot(x0,x3,lw=3,c='r',label='$2\\nu$ probability')
#ax1.plot(xx0,xx1,marker='o',c='b
major_ticksx = np.arange(0, 30,1)
minor_ticksx = np.arange(0, 30,0.5)
#major_ticksy = np.arange(0, 0.06,0.01)
#minor_ticksy = np.arange(0, 0.06, 0.005)
ax1.set_xticks(major_ticksx)
ax1.set_xticks(minor_ticksx, minor=True)
#ax1.set_yticks(major_ticksy)
#ax1.set_yticks(minor_ticksy, minor=True)
ax1.tick_params(axis='both', which='major', labelsize=20,width=1,length=6)
ax1.tick_params(axis='both', which='minor',width=1,length=3)
ax1.set_xlim([12,20])
#ax1.set_ylim([0,0.05])
leg = ax1.legend(handlelength=2)
leg = ax1.legend()
#leg = ax2.legend()
plt.tight_layout()
fig.savefig('current-Q-v3.pdf', bbox_inches='tight')
plt.show()

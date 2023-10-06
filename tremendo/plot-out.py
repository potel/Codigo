#!/usr/bin/env python3
import numpy as np
import math
import matplotlib.pyplot as plt
import sys
from matplotlib import rc,rcParams
import random as ran
rc('text',usetex=True)
rc('font',**{'family':'serif','serif':['Computer Modern'],'size':13})
def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False
color=['k', 'r', 'b', 'g', 'm','y','c','burlywood','bisque','lime','orange',
       'tomato','deeppink','indigo']
numargs=len(sys.argv)
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.set_yscale('log')
ax1.set_xlabel('$\\theta_{cm}$',fontsize=20)
ax1.set_ylabel('$d\sigma/d\Omega$',fontsize=20)
count=0
with open(sys.argv[1]) as f1:
    lines = f1.readlines()
    titlefig=f"{lines[3]}"
    print(titlefig)
    for i in range(len(lines)):
        if(lines[i].strip()):
            if (is_number(lines[i].split()[0])):
                count=count+1
numlines=count
print('Number of lines: ',numlines)
x0=[0 for j in range(count)]
x1=[0 for j in range(count)]
count=0
with open(sys.argv[1]) as f1:
    lines = f1.readlines()
    for i in range(len(lines)):
        if(lines[i].strip()):
            if (is_number(lines[i].split()[0])):
                x0[count]=(float(lines[i].split()[0]))
                x1[count]=(float(lines[i].split()[1]))
                count=count+1
ax1.plot(x0,x1,lw=2,c='k',label='first set')
for j in range(1,numargs-1):
    count=0
    with open(sys.argv[j+1]) as f1:
        lines = f1.readlines()
        for i in range(len(lines)):
            if(lines[i].strip()):
                if (is_number(lines[i].split()[0])):
                    x0[count]=(float(lines[i].split()[0]))
                    x1[count]=(float(lines[i].split()[1]))
                    count=count+1
    st=f'set {j}'
    ii=ran.randint(0,len(color)-1)
    ax1.plot(x0,x1,lw=2,c=color[ii],label=st)


#ax1.plot(xx0,xx1,marker='o',c='b
major_ticksx = np.arange(0, 190, 20)
minor_ticksx = np.arange(0, 190, 10)
#major_ticksy = np.arange(0, 5, 0.5)
#minor_ticksy = np.arange(0, 10, 0.1)
ax1.set_xticks(major_ticksx)
ax1.set_xticks(minor_ticksx, minor=True)
#ax1.set_yticks(major_ticksy)
#ax1.set_yticks(minor_ticksy, minor=True)
ax1.tick_params(axis='both', which='major', labelsize=20,width=1,length=6)
ax1.tick_params(axis='both', which='minor',width=1,length=3)
ax1.set_xlim([0,70])
#ax1.set_ylim([1.3,2.6])
leg = ax1.legend(handlelength=2)
leg = ax1.legend()
#leg = ax2.legend()
plt.tight_layout()
fig.savefig('fig-out-c3.pdf', bbox_inches='tight')
plt.show()

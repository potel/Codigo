import numpy as np
import matplotlib.pyplot as plt
import sys
csfont = {'fontname':'Times'}
hfont = {'fontname':'Helvetica'}
#print sys.argv[0]
#print sys.argv[1]

with open(sys.argv[1]) as f:
    lines = f.readlines()
    x = [line.split()[int(sys.argv[2])] for line in lines]
    y = [line.split()[int(sys.argv[3])] for line in lines]
fig = plt.figure()

ax1 = fig.add_subplot(111)
#ax2 = fig.add_subplot(122)

ax1.set_title(r"$^{128}$Te ($^3$He,$n$)$^{130}$Xe, $E_{lab}=24 MeV$",**csfont)    
ax1.set_xlabel(r"$\theta_{CM}$")
ax1.set_ylabel('$d\sigma/d\Omega (\mu b/sr)$')

ax1.semilogy(x,y, c='b',lw=3, label='$0^+$ ground state')
#ax2.plot(x,y2, c='b', label='hey2')
leg = ax1.legend()
#leg = ax2.legend()

plt.show()

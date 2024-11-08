from pylab import *
from thermo.gpumd.data import load_kappa, load_shc
from thermo.gpumd.calc import running_ave, hnemd_spectral_kappa
import pandas as pd
from glob import glob


aw = 2
fs = 16
font = {'size'   : fs}
matplotlib.rc('font', **font)
matplotlib.rc('axes' , linewidth=aw)

def set_fig_properties(ax_list):
    tl = 8
    tw = 2
    tlm = 4
    
    for ax in ax_list:
        ax.tick_params(which='major', length=tl, width=tw)
        ax.tick_params(which='minor', length=tlm, width=tw)
        ax.tick_params(axis='both', direction='in', right=True, top=True)


dt = 0.001 #ps
Ns = 1000 #Sample interval
thermo = loadtxt("compute_Energy_Temp.out")
# print(thermo)
jp = loadtxt("compute_Heatflux.out")
# print(jp)
BLOCK_LENGTH = 426
start = 250
Ein = thermo[start:, 1]
Eout = thermo[start:, 2]
# print(Ein)
# print(Eout)
Etol=(Eout - Ein) / 2 / 1000 #in units of Kev
# print(Etol)
Etol = Etol - Etol[0]
# print(Etol)
t = dt * np.arange(1, len(Etol) + 1) * Ns / 1000 #unit in ns
# print(t)
jpy = jp[start:, 2] - jp[start:, 5]

jpy = jpy / BLOCK_LENGTH / 10 * 1000 #in units of ev/ns
# print(jpy)
accum_jpy = cumsum(jpy) * 0.001 / 1000 #in units of kev
# print(cumsum(jpy))
# print(t)
# print(accum_jpy)
# print(Etol)
with open('compare_modified.txt', 'w') as fl:
    for i in range(0, len(t)):
        fl.writelines('{:.3f}'.format(t[i]) + ' ')
        fl.writelines('{:.6f}'.format(np.average(abs(accum_jpy[i]))) + ' ')
        fl.writelines('{:.6f}'.format(np.average(abs(Etol[i]))) + '\n')

figure(figsize=(6,5))
# subplot(1,2,1)
set_fig_properties([gca()])
plot(t, abs(accum_jpy), 'C3', linewidth=3, label = 'from atoms')
plot(t, abs(Etol), 'C0', linewidth=3, label = 'from thermostats')
# gca().set_xticks(group_idx)
xlabel('Time (ns)')
ylabel('Accumulated heat (keV)')
legend()
savefig('heatflux.png')
show()
close()


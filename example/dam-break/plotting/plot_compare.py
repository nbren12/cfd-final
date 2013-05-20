from pylab import *
from mpl_toolkits.mplot3d import Axes3D

run = 'dam_break1d_n2048_superbee'
run_hr = 'dam_break1d_n64_superbee'
run_mm = 'dam_break1d_n64_minmod'
run_o1 = 'dam_break1d_n64_o1'

frame= 7
hr = __import__(run)
hr.cont.read_frame(frame)

xx = hr.cont.state.grid.x.centers
ex = hr.cont.state.q[0,:,0]

fig = figure()
ax = fig.add_subplot(111)


ll = [run_hr,run_mm,run_o1]
data = map(__import__,ll)
map(lambda x : x.cont.read_frame(frame),data)
x = data[0].cont.state.grid.x.centers

hhr, hmm,h01 = map(lambda x: x.cont.state.q[0,:,0],data)


ax.plot(xx,ex,'k')
ax.plot(x,hhr,'ko')
ax.plot(x,hmm,'g-')
ax.plot(x,h01,'b')
ax.set_xlim([-7,7])
ax.set_ylim([.95,3.05])
ax.set_xlabel('x')
ax.set_ylabel('h')


plt.legend(['n = 2048 Superbee','n = 64 Superbee','n = 64 MinMod','n = 64 Godunov'])
fig.savefig('t%.1f_compare.eps'%data[0].cont.time,bbox_inches=None)

show()

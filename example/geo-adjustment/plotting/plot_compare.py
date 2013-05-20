from pylab import *
from mpl_toolkits.mplot3d import Axes3D

run = 'geo_break1d_n1024_hr'
run_hr = 'geo_break1d_n128_hr'
run_o1 = 'geo_break1d_n128_o1'

frame= 6
hr = __import__(run)
hr.cont.read_frame(frame)

xx = hr.cont.state.grid.x.centers/10
ex = hr.cont.state.q[0,:,0]

fig = figure()
ax = fig.add_subplot(111)


ll = [run_hr,run_o1]
data = map(__import__,ll)
map(lambda x : x.cont.read_frame(frame),data)
x = data[0].cont.state.grid.x.centers/10

hhr, h01 = map(lambda x: x.cont.state.q[0,:,0],data)


ax.plot(xx,ex,'k')
ax.plot(x,hhr,'ko')
# ax.plot(x,hmm,'g-')
ax.plot(x,h01,'b')
# ax.set_xlim([-7,7])
# ax.set_ylim([.95,3.05])
ax.set_xlabel('x')
ax.set_ylabel('h')


plt.legend(['n = 1024 HR','n = 128 HR','n = 128 Godunov'])
fig.savefig('t%.1f_compare.eps'%data[0].cont.time,bbox_inches=None)

show()

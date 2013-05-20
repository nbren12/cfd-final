from pylab import *
from mpl_toolkits.mplot3d import Axes3D


run = 'geo_break1d_n1024_hr'

data = __import__(run)

fig, axs = plt.subplots(nrows=3,ncols=2,figsize=(4,7))

for i in range (3):
    for j in range(2):
        data.cont.read_frame(i*2+j)
        x = data.cont.state.grid.x.centers
        h = data.cont.state.q[0,:,0]
        ax = axs[i,j]
        ax.plot(x,h)
        # ax.set_ylim([.9,3.05])
        # ax.set_xlim([-7,7])
        ax.set_title('t=%.1f'%data.cont.time)
        if i != 3:
            ax.set_xticks([])
        else:
            ax.set_xlabel('x')
        if j!= 0:
            ax.set_yticks([])
        else:
            ax.set_ylabel('h')


show()

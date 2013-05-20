from pylab import *
from mpl_toolkits.mplot3d import Axes3D


run = 'hr1d_n4000'

data = __import__(run)

fig, axs = plt.subplots(nrows=4,ncols=2,figsize=(4,7))

for i in range (4):
    for j in range(2):
        data.cont.read_frame(i*2+j)
        x = data.cont.state.grid.x.centers
        h = data.cont.state.q[0,:,0]
        ax = axs[i,j]
        ax.plot(x,h)
        ax.set_ylim([.5,3.05])
        ax.set_xlim([0,7])
        ax.set_title('t=%.1f'%data.cont.time)
        if i != 3:
            ax.set_xticks([])
        else:
            ax.set_xlabel('r')
        if j!= 0:
            ax.set_yticks([])
        else:
            ax.set_ylabel('h')
fig.savefig('rad-plot-evol.eps',bbox_inches=None)

# run = 'n500_2dhr'
# data = __import__(run)
#
# fig1 = plt.figure()#   subplots(nrows=4,ncols=2,figsize=(4,7))
# ax = fig1.add_subplot(111,projection='3d')
# data.cont.read_frame(4)
# data.cont.surf_plot(ax=ax)
# fig1.savefig('t%.1f_2dhr.eps'%data.cont.time,bbox_inches=None)

# for i in range (4):
#     for j in range(2):
#         ax = fig1.add_subplot(4,2,i*2+j,projection='3d')
#         data.cont.read_frame(i*2+j)
#         x,y = data.cont.state.grid.c_centers
#         h = data.cont.state.q[0,:,:]
#         data.cont.surf_plot(ax=ax)
#         # ax.set_title('t=%.1f'%data.cont.time)
#         # if i != 3:
#         #     ax.set_xticks([])
#         # else:
#         #     ax.set_xlabel('r')
#         # if j!= 0:
#         #     ax.set_yticks([])
#         # else:
#         #     ax.set_ylabel('h')
#
show()

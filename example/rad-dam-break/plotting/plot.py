import numpy as np
from clawpack import pyclaw
import os
from swe2d import Controller, ControllerSW2D,advance_sw
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


run = 'n100_hr_nocfix'
# run = 'n100_2dhr'

mod = __import__(run)

cont = mod.cont
cont.read_frame(2)
h = cont.state.q[0,:,:]

# cont.read_frame(0)
# fig = plt.figure(figsize=(12,4))
# ax = fig.add_subplot(121,projection='3d')
# cont.surf_plot(ax=ax)
#
# cont.read_frame(7)
fig = plt.figure()

ax = fig.add_subplot(111,projection='3d')
cont.surf_plot(ax)
fig.savefig('%s_t%f_3d.eps'%(run,cont.time),bbox_inches=None)



fig1 = plt.figure()
ax = fig1.add_subplot(111)
cont.cont_plot(ax=ax)
ax.set_aspect('equal')
# fig1.savefig('%s_t%f.eps'%(run,cont.time),bbox_inches=None)
plt.show()


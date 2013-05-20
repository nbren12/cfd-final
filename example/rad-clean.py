"""
Shallow water Equations in One Dimension with no bottom topography are:

h_t + (hu)_x + (hv)_y = 0
(hu)_t + (hu^2 + 1/2 g h^2 )_x + (huv)_y = 0
(hv)_t + (huv)_x + (h v^2 + 1/2 g h^2)_y = 0

hu^2 = q2^2/q1 + 1/2 g q2^2

This is obviously symmetric in x and y.
Here, I used periodic boundary conditions. Might want to generalize the code later.
"""
import numpy as np
from clawpack import pyclaw
import os
from swe2d import Controller, ControllerSW2D,advance_sw
from matplotlib import pyplot as plt
from scipy.interpolate import RectBivariateSpline

suffix = "hr_nocfix"

R = 3
T  = 1.0
g = 9.812           # Gravity
H = 2               # Average Depth
eta = 1             # Height Deviation
c = np.sqrt(H*g)    # Speed of Gravity waves

# Initialize Grid
n = 1000
nx =n
ny =n
folder = "n%d_%s"%( nx,suffix )


#   First Dimension is x, Second Dimension is y
x = pyclaw.Dimension('x',-10.0,10.0,nx)
y = pyclaw.Dimension('y',-10.0,10.0,ny)
domain = pyclaw.Domain((x,y))

# Initialized state and Problem Parameters
state  = pyclaw.State(domain,3)
state.problem_data ={  'g':g , 'efix':True,'hr':True,'bcs':0,'cfix':0}
s_opts = {'f':0.1}


# Figure out the time step
dx,dy = state.grid.delta
dt = min(dx,dy)/c/1.5
nt = int(T/dt)


print("Beginning Time Stepping")
cont = ControllerSW2D(state,advance_sw,dt=dt,prefix=folder)
cont.run(T)

# from mpl_toolkits.mplot3d import Axes3D
# h = cont.state.q[0,:,:]
#
# cont.read_frame(0)
# fig = plt.figure(figsize=(12,4))
# ax = fig.add_subplot(121,projection='3d')
# cont.surf_plot(ax=ax)
#
# cont.read_frame(5)
# ax = fig.add_subplot(122,projection='3d')
# cont.surf_plot(ax)
#
#
# fig1 = plt.figure()
# ax = fig1.add_subplot(111)
# cont.cont_plot(ax=ax)
# plt.show()
#

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
from swe2d import Controller, advance_sw,advance_coriolis
from matplotlib import pyplot as plt



g = .1          # Gravity
H = 10               # Average Depth
eta = H*.01            # Height Deviatio
c = np.sqrt(H*g)    # Speed of Gravity waves
T  = 1.0/c/1.5
f  = 1.0/c *5

suffix = "o2"

# n = 7
n = 7

nx = 2**n
ny = 2


#   First Dimension is x, Second Dimension is y
x = pyclaw.Dimension('x',-1.0,1.0,nx)
y = pyclaw.Dimension('y',-1.0,1.0,ny)

domain = pyclaw.Domain((x,y))
state  = pyclaw.State(domain,3)
state.problem_data ={  'g':g , 'efix':True,'hr':True,'bcs':2,'cfix':0}
s_opts = {'f':f}

dx,dy = state.grid.delta
dt = min(dx,dy)/c/3.0
nt = int(T/dt)


# Initial Data for 1d Dam Break
h0 = np.ones((nx,ny))*H#
h0 -=np.sign(x.centers)[:,None]*eta
u0 = np.zeros(domain.grid.num_cells)
v0 = np.zeros(domain.grid.num_cells)


state.q[0,:,:] = h0
state.q[1,:,:] = (u0*h0)
state.q[2,:,:] = (v0*h0)

cont = Controller(state,advance_sw,dt=dt,srcsolver=(advance_coriolis,s_opts) ,prefix="geo_break1d_n%d_%s"%(nx,suffix))

cont.run(T)

xx,yy,z = cont.get_plot_args()
plt.plot(xx[:,0],z[:,0])
# plt.contourf(*cont.get_plot_args())
plt.show()


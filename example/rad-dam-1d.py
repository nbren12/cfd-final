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
from swe2d import Controller, advance_sw,advance_metric_terms
from matplotlib import pyplot as plt

R = 3
T  = 1.0
g = 9.812           # Gravity
H = 2               # Average Depth
eta = 1             # Height Deviation
c = np.sqrt(H*g)    # Speed of Gravity waves

# Initialize Grid
nx = 4000
ny = 2

x = pyclaw.Dimension('x',0.0,15.0,nx)
y = pyclaw.Dimension('y',0.0,10.0,ny)
domain = pyclaw.Domain((x,y))

# Initialized state and Problem Parameters
state  = pyclaw.State(domain,3)
state.problem_data ={  'g':g , 'efix':True,'hr':True,'bcs':3,'cfix':0}
s_opts = {'r':state.grid.c_centers[0]}

# Initial Data for 1d Dam Break
h0 = np.ones((nx,ny))*H#
h0 -=np.sign(x.centers-R)[:,None]*eta
u0 = np.zeros(domain.grid.num_cells)
v0 = np.zeros(domain.grid.num_cells)

# Figure out the time step
dx,dy = state.grid.delta
dt = min(dx,dy)/c/1.5
nt = int(T/dt)


state.q[0,:,:] = h0
state.q[1,:,:] = (u0*h0)
state.q[2,:,:] = (v0*h0)

cont = Controller(state,advance_sw,srcsolver=(advance_metric_terms,s_opts),dt=dt,prefix="hr1d_n4000_ss")
cont.run(T)

xx,yy,z = cont.get_plot_args()
plt.plot(xx[:,0],z[:,0])
# plt.contourf(*cont.get_plot_args())
plt.show()


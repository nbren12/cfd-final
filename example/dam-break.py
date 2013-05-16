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
from swe2d import Controller, advance_sw
from matplotlib import pyplot as plt



g = 9.812
H = 1
nx = 100
ny = 2


#   First Dimension is x, Second Dimension is y
x = pyclaw.Dimension('x',-10.0,10.0,nx)
y = pyclaw.Dimension('y',-10.0,10.0,ny)

domain = pyclaw.Domain((x,y))
state  = pyclaw.State(domain,3)
state.problem_data ={  'g':g , 'efix':True,'hr':True,'bcs':2}
s_opts = {'f':0.1}

dx,dy = state.grid.delta
dt = dx / 10
T  = 0.1
nt = int(T/dt)


# Initial Data for 2d Radial Dam Break
h0 = np.ones((nx,ny)) #
h0 = (2-np.sign(x.centers))[:,None]
u0 = np.zeros(domain.grid.num_cells)
v0 = np.zeros(domain.grid.num_cells)


state.q[0,:,:] = h0
state.q[1,:,:] = (u0*h0)
state.q[2,:,:] = (v0*h0)

cont = Controller(state,advance_sw)
cont.run(T)

xx,yy,z = cont.get_plot_args()
plt.plot(xx,z)
# plt.contourf(*cont.get_plot_args())
plt.show()


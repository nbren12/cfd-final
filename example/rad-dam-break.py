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
ny = 100


#   First Dimension is x, Second Dimension is y
x = pyclaw.Dimension('x',-10.0,10.0,nx)
y = pyclaw.Dimension('y',-10.0,10.0,ny)

domain = pyclaw.Domain((x,y))
state  = pyclaw.State(domain,3)
state.problem_data ={  'g':g , 'efix':True,'hr':True,'bcs':0}
s_opts = {'f':0.1}

dx,dy = state.grid.delta
dt = dx / 10
T  = .1
nt = int(T/dt)


# Initial Data for 2d Radial Dam Break
h0 = np.ones((nx,ny)) #
h0 += 2*(reduce(lambda x,y:np.sqrt(x+y),map(lambda x: x**2,domain.grid.c_centers))< 4)
# h0 += np.exp(-reduce(lambda x,y:np.sqrt(x+y),map(lambda x: x**2,domain.grid.c_centers)))
u0 = np.zeros(domain.grid.num_cells)
v0 = np.zeros(domain.grid.num_cells)


state.q[0,:,:] = h0
state.q[1,:,:] = (u0*h0)
state.q[2,:,:] = (v0*h0)

cont = Controller(state,advance_sw)
cont.run(T)

from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax = fig.add_subplot(111,projection='3d')
ax.plot_surface(*cont.get_plot_args(),
        rstride=1, cstride=1,
        color='green',linewidth=.1,shade=True)

# plt.contourf(*cont.get_plot_args())
plt.show()


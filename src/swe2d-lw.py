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
from ipdb import set_trace as st
from math import sqrt
from clawpack import pyclaw
from matplotlib import pyplot as plt
from hr_solver2d import advance_1d

g = 9.812
nx = 100
ny = 100
ng = 2

#   First Dimension is x, Second Dimension is y
x = pyclaw.Dimension('x',-10.0,10.0,nx)
y = pyclaw.Dimension('y',-10.0,10.0,nx)

dx,dy = state.grid.delta
dt = dx / 10
T  = 1.0
nt = int(T/dt)

domain = pyclaw.Domain((x,y))
state  = pyclaw.State(domain,3)

# Initial Data for 1D Riemann Problem
h0 = (2-np.sign(x.centers))
u0 = np.zeros(x.num_cells)
v0 = np.zeros(x.num_cells)

state.q[0,:,:] = h0[:,None]
state.q[1,:,:] = (u0*h0)[:,None]
state.q[2,:,:] = (v0*h0)[:,None]

for i in xrange(1):
    # qbc = advance_1d(state.q,g,dt,dx)
    pass

# plt.plot(state.grid.c_centers[0],state.q[0,:])
# plt.show()

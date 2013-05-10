"""
Shallow water Equations in One Dimension with no bottom topography are:

h_t + (hu)_x = 0
(hu)_t + (hu^2 + 1/2 g h^2 )_x = 0

hu^2 = q2^2/q1 + 1/2 g q2^2

Here, I used periodic boundary conditions. Might want to generalize the code later.
"""
import numpy as np
from ipdb import set_trace as st
from math import sqrt
from clawpack import pyclaw
from matplotlib import pyplot as plt
from hr_solver1d import advance_1d

nx = 500
ng = 2
x = pyclaw.Dimension('x',-10.0,10.0,nx)
domain = pyclaw.Domain(x)
state = pyclaw.State(domain,2)

# Initial Data for 1D Riemann Problem
h0 = (2-np.sign(x.centers))
u0 = np.zeros(x.num_cells)

state.q[0,:]=h0
state.q[1,:]=u0*h0

state.problem_data['grav'] = 9.81
state.problem_data['efix'] = False

g = state.problem_data['grav']
dx = state.grid.delta[0]
dt = dx / 10
T  = 1.0
nt = int(T/dt)
for i in xrange(nt):
    qbc = advance_1d(state.q,g,dt,dx)
    # pass

plt.plot(state.grid.c_centers[0],state.q[0,:])
plt.show()

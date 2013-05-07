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
from calc_flux_roe1d_lw import calc_flucts,advance_1d
from calc_flux_roe1d import calc_fluxes
nx = 200
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


qbc = np.zeros((2,ng*2 + nx ))
fluxes = np.zeros((2,nx+1))
g = state.problem_data['grav']
dx = state.grid.delta[0]
dt = dx / 10
T  = .5
nt = int(T/dt)
for i in xrange(nt):
# for i in xrange(1):
    print i

    # fill ghost cell
    qbc[:,ng:-ng] = state.q
    qbc[:,:ng] = state.q[:,-ng:]
    qbc[:,-ng:] = state.q[:,:ng]

    # Calculate Fluxes at every cell
    alpha,s,R  = calc_flucts(qbc,nx,ng,g,dt,dx)
    F,F_tilde = advance_1d(qbc,nx,ng,g,dt,dx)
    state.q = state.q + dt/dx * (F +F_tilde[:,1:]-F_tilde[:,:-1])
    # Advance Solution
    # state.q =state.q -(dt/dx)*(fluxes[:,1:]-fluxes[:,:-1])

plt.plot(state.grid.c_centers[0],state.q[0,:])
plt.show()

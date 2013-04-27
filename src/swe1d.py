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
nx = 512
ng = 1
x = pyclaw.Dimension('x',-10.0,10.0,nx)
domain = pyclaw.Domain(x)
state = pyclaw.State(domain,2)

# Initial Data for 1D Riemann Problem
num_ghost = 1
h0 = (2-np.sign(x.centers))
u0 = np.zeros(x.num_cells)

state.q[0,:]=h0
state.q[1,:]=u0*h0

state.problem_data['grav'] = 9.81
state.problem_data['efix'] = False

def calc_flux(ql,qr,g):

    ul = ql[1]/ql[0]
    ur = qr[1]/qr[0]
    hl = ql[0]
    hr = qr[0]

    hl2 = sqrt(hl)
    hr2 = sqrt(hr)

    fql = np.array([hl*ul,hl*ul**2+g*hl**2/2])
    fqr = np.array([hr*ur,hr*ur**2+g*hr**2/2])

    u_hat = ( hl2*ul + hr2*ur )/(hl2 + hr2)
    h_bar = (hl+hr)/2
    c_hat = sqrt(g*h_bar)



    A_hat = np.matrix([[0,1],[-u_hat**2 + g * h_bar,2*u_hat]])

    lam_hat = np.matrix(np.diag([u_hat-c_hat,u_hat+c_hat]))
    R = np.matrix([[1,1],[u_hat-c_hat,u_hat+c_hat]])
    L = np.matrix([[u_hat+c_hat,-1],[-u_hat+c_hat,1]]) / 2.0 / c_hat

    A_abshat = R * np.abs(lam_hat) * L
    flux = (fql+fqr)/2.0  - np.dot(A_abshat,qr-ql)/2.0

    return flux

qbc = np.zeros((2,ng*2 + nx ))
fluxes = np.zeros((2,nx+1))
g = state.problem_data['grav']
dx = state.grid.delta[0]
dt = dx / 10
T  = 1.0
nt = int(T/dt)
for i in xrange(nt):

# fill ghost cell
    qbc[:,ng:-ng] = state.q
    qbc[:,0] = state.q[:,-1]
    qbc[:,-1]= state.q[:,0]

# Calculate Fluxes at every cell
    for j in xrange(nx+ng):
        fluxes[:,j]  = calc_flux(qbc[:,j],qbc[:,j+1],g)
# Advance Solution
    state.q =state.q -(dt/dx)*(fluxes[:,ng:]-fluxes[:,:-ng])

plt.plot(state.grid.c_centers[0],state.q[0,:])
plt.show()

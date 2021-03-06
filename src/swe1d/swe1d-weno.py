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
from weno_solver1d import advance_1d
from weno_f2py import weno
from pylab import *

nx = 50
ng = 3
x = pyclaw.Dimension('x',-10.0,10.0,nx)

domain = pyclaw.Domain(x)
state = pyclaw.State(domain,2)
q0 = np.load('q.npy')
# Initial Data for 1D Riemann Problem
h0 = (2-np.sign(x.centers))
u0 = np.zeros(x.num_cells)

state.q[0,:]=h0
state.q[1,:]=u0*h0
# state.q = q0

state.problem_data['grav'] = 9.81
state.problem_data['efix'] = False

g = state.problem_data['grav']
dx = state.grid.delta[0]
dt = dx / 10
T  = 1.0
nt = int(T/dt)
ion()
figure()
for i in xrange(nt):
    if i>0 : qm1= state.q
    state.q,q_pm,qbc,qb =advance_1d(state.q,g,dt,dx)
    clf()
    plt.plot(x.centers, qb[0,:])
    axis([-10,10,1,3])
    plt.pause(.1)
    ql = q_pm[:,:,0]
    qr = q_pm[:,:,1]

# plt.plot(x.centers,qm1[0,:],'k')
# plt.plot(x.edges,ql[0,:],'o')
# plt.plot(x.edges,qr[0,:],'go')
# plt.plot(qb[0,:])
# plt.plot(state.q[0,:])
# plt.show()

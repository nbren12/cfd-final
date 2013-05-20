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
from adv_weno import advance_1d
from pylab import *

nx = 100
ng = 3
x = pyclaw.Dimension('x',-1/2.0,1/2.0,nx)

domain = pyclaw.Domain(x)
state = pyclaw.State(domain,1)

# Initial Data for 1D Riemann Problem
u_ex = lambda T : cos(pi*(x.centers-T))**(10) * (1+sign(x.centers-T))/2
u0 = zeros((1,nx))
u0[0,:] = u_ex(0)

q= u0.copy()

dx = state.grid.delta[0]
dt = dx / 2
T  = 1.0
nt = int(T/dt)
figure()
for i in xrange(nt):

    L = advance_1d(q,dt,dx)
    qs = q + dt*L
    L = advance_1d(qs,dt,dx)
    qs = (3.0/4.0) *q + (1.0/4.0) *qs + (1.0/4.0)*dt*L
    L = advance_1d(qs,dt,dx)
    qs = (1.0/3.0)*q  + (2.0/3.0) *qs + (2.0/3.0) *dt*L
    q= qs

    # L = advance_1d(q,dt,dx)
    # q = q + dt*L

clf()
# plot(x.centers,u_ex(T))
plot(x.centers,u0[0,:])
plot(x.centers, q.T,'o')
# plot(x.edges,ql.T,'ro')
# plot(x.edges,qr.T,'bo')
# plot(x.edges,q_pm[0,:,0],'bo')
# plot(x.edges,q_pm[0,:,1],'ro')
axis([-.5,.5,0,1])
# plt.pause(.1)
show()


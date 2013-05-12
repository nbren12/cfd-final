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
from rp_sw2d_roe_f2py import rp_sw2d_roe as rp

g = 9.812
nx = 100
ny = 100
ng = 2

nframes = 10

#   First Dimension is x, Second Dimension is y
x = pyclaw.Dimension('x',-10.0,10.0,nx)
y = pyclaw.Dimension('y',-10.0,10.0,nx)

domain = pyclaw.Domain((x,y))
state  = pyclaw.State(domain,3)

dx,dy = state.grid.delta
dt = dx / 10
T  = 1.0
nt = int(T/dt)


# Initial Data for 1D Riemann Problem
# h0 = (2-np.sign(x.centers))[:,None] # Dam Break
h0 = np.ones((nx,ny))
# h0[nx/2-nx/3:nx/2+nx/3,ny/2-ny/3:ny/2+ny/3] = 3
# h0[40:60,40:60] = 3
# h0 += np.exp(-reduce(lambda x,y:np.sqrt(x+y),map(lambda x: x**2,domain.grid.c_centers)))
h0 += (reduce(lambda x,y:np.sqrt(x+y),map(lambda x: x**2,domain.grid.c_centers))< 1)
u0 = np.zeros(domain.grid.num_cells)
v0 = np.zeros(domain.grid.num_cells)


xx,yy = state.grid.c_centers
xe,ye = state.grid.c_edges
# state.q[0,:,:] = h0[:,None]
state.q[0,:,:] = h0
state.q[1,:,:] = (u0*h0)
state.q[2,:,:] = (v0*h0)


fig, axl = plt.subplots(ncols=nframes/2,nrows=2,figsize=(12,5))
axl = axl.ravel()

dframe = nt/nframes
frame= 0

for i in xrange(nt):
    if i%dframe == 0:
        axl[frame].contourf(xx,yy,state.q[0,:,:])
        axl[frame].set_aspect('equal')
        axl[frame].set_title('T = %.2f'%(i*dt))
        # axl[i].
        frame+=1
    qm = state.q
    advance_1d(state.q,g,dt,dx,dy,2)


plt.show()

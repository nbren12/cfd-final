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
from hr_solver2d import advance_sw, advance_coriolis

g = .1
ymult = 3
xmult = 3
ny = ymult*32
nx = xmult*32
ng = 2

H = 10.0
c = np.sqrt(g*H)
L = 100.0
T  = L/c * .2
# T = .5
eta = H*.3
f = L/c * 0

xL = xmult*L
yL = ymult*L
nframes = 10

#   First Dimension is x, Second Dimension is y
x = pyclaw.Dimension('x',-xmult*L,xmult*L,nx)
y = pyclaw.Dimension('y',-ymult*L,ymult*L,ny)

domain = pyclaw.Domain((x,y))
state  = pyclaw.State(domain,3)

dx,dy = state.grid.delta
dt = dx / 10
nt = int(T/dt)


# Initial Data for 1D Riemann Problem
h0 = H   -eta*np.sign(x.centers)[:,None] # Dam Break
# h0 = np.ones((nx,ny))
# h0[nx/2-nx/3:nx/2+nx/3,ny/2-ny/3:ny/2+ny/3] = 3
# h0[40:60,40:60] = 3
# h0 += np.exp(-reduce(lambda x,y:np.sqrt(x+y),map(lambda x: x**2,domain.grid.c_centers)))
# h0 += (reduce(lambda x,y:np.sqrt(x+y),map(lambda x: x**2,domain.grid.c_centers))< 3)
u0 = np.zeros(domain.grid.num_cells)
u0[:,ny/2-ny/6:ny/2+ ny/6] = c * eta
v0 = np.zeros(domain.grid.num_cells)


xx,yy = state.grid.c_centers
xe,ye = state.grid.c_edges
state.q[0,:,:] = h0
state.q[1,:,:] = (u0*h0)
state.q[2,:,:] = (v0*h0)



fig, axl = plt.subplots(ncols=nframes/2,nrows=2,figsize=(12,5))
fig1, axl1 = plt.subplots(ncols=nframes/2,nrows=2,figsize=(12,5))
axl = axl.ravel()
axl1 = axl1.ravel()

pltfun = lambda q: q[1,:,:]/q[0,:,:]
# pltfun = lambda q: q[0,:,:]
mmax =  c *eta*1.3
mmin =  -c * eta*1.3
dframe = nt/nframes
frame= 0

for i in xrange(nt):
    if i%dframe == 0 and frame <10:
        # axl[frame].contourf(xx,yy,pltfun(state.q),np.linspace(mmin,mmax,12), vmin=mmin,vmax=mmax)
        axl[frame].pcolor(xx,yy,pltfun(state.q))
        axl[frame].set_aspect('equal')
        axl[frame].set_xlim([-L,L])

        axl1[frame].plot(x.centers,pltfun(state.q)[:,ny/2],'r')
        axl1[frame].set_xlim([-L,L])
        axl1[frame].set_ylim([mmin,mmax])
        axl1[frame].set_title('T = %.2f'%(i*dt))
        if frame <5:
            axl[frame].get_xaxis().set_visible(False)
            axl1[frame].get_xaxis().set_visible(False)
        if frame not in [0,5]:
            axl[frame].get_yaxis().set_visible(False)
            axl1[frame].get_yaxis().set_visible(False)
        frame+=1
    qm = state.q
    print "%d of %d"%(i,nt)
    advance_sw(state.q,g,dt,dx,dy)
    advance_coriolis(state.q,f,dt)

plt.show()

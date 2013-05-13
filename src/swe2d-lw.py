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
from hr_solver2d import advance_sw

import os
from SaveVideo import *
from PIL import Image

g = .1
H = 100
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
T  = 25
nt = int(T/dt)


# Initial Data for 1D Riemann Problem
# h0 = (2-np.sign(x.centers))[:,None] # Dam Break
h0 = np.ones((nx,ny))*H #   + np.sign(x.centers)*H*.1
# h0[nx/2-nx/3:nx/2+nx/3,ny/2-ny/3:ny/2+ny/3] = 3
# h0[40:60,40:60] = 3
# h0 += 2*np.exp(-reduce(lambda x,y:np.sqrt(x+y),map(lambda x: x**2,domain.grid.c_centers)))
# h0 += 2*(reduce(lambda x,y:np.sqrt(x+y),map(lambda x: x**2,domain.grid.c_centers))< 4)
# u0 = np.zeros(domain.grid.num_cells)
u0 = -np.ones(domain.grid.num_cells)*np.sqrt(g*H)*2
u0[:,ny/2:ny/2+ny/6] = np.sqrt(g*H)*2
u0[:,ny/2:] = np.sqrt(g*H)*2
v0 = np.zeros(domain.grid.num_cells)


xx,yy = state.grid.c_centers
xe,ye = state.grid.c_edges
state.q[0,:,:] = h0
state.q[1,:,:] = (u0*h0)
state.q[2,:,:] = (v0*h0)


# fig, axl = plt.subplots(ncols=nframes/2,nrows=2,figsize=(12,5))
# fig1, axl1 = plt.subplots(ncols=nframes/2,nrows=2,figsize=(12,5))
# axl = axl.ravel()
# axl1 = axl1.ravel()

vs = VideoSink((nx,ny),filename='bone')
from matplotlib import cm
cmm = cm.ScalarMappable(cmap=cm.jet)

mmin = 0
mmax = 2
dframe = nt/nframes
frame= 0
fig = plt.figure()
j = 0
for i in xrange(nt):
    if i%5 ==0:
        u = state.q[1,:,:]/state.q[0,:,:]
        v = state.q[2,:,:]/state.q[0,:,:]
        cmm.set_array(u)
        cmm.autoscale()

        # w = (v[1:,1:] -v[:-1,1:])/dy -(u[1:,1:]-u[1:,:-1])/dx
        writeme = cmm.to_rgba(u,bytes=True)
        pil = Image.fromarray(writeme)
        vs.run(pil)

    # if i%dframe == 0 and frame <10:
    #     axl[frame].contourf(xx,yy,state.q[0,:,:],np.linspace(mmin,mmax,12), vmin=mmin,vmax=mmax)
    #     axl[frame].set_aspect('equal')

    #     axl1[frame].plot(x.centers,state.q[0,:,ny/2],'r')
    #     axl1[frame].set_xlim([-10,10])
    #     axl1[frame].set_ylim([mmin,mmax])
    #     axl1[frame].set_title('T = %.2f'%(i*dt))
    #     if frame <5:
    #         axl[frame].get_xaxis().set_visible(False)
    #         axl1[frame].get_xaxis().set_visible(False)
    #     if frame not in [0,5]:
    #         axl[frame].get_yaxis().set_visible(False)
    #         axl1[frame].get_yaxis().set_visible(False)
    #     # axl[i].
    #     frame+=1
    # qm = state.q
    advance_sw(state.q,g,dt,dx,dy,efix=True,hr=True)
vs.close()
# plt.show()

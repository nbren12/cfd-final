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
from scipy.interpolate import RectBivariateSpline

R = 3
T  = .2
g = 9.812           # Gravity
H = 2               # Average Depth
eta = 1             # Height Deviation
c = np.sqrt(H*g)    # Speed of Gravity waves

# Initialize Grid
nx = 50
ny = 50


#   First Dimension is x, Second Dimension is y
x = pyclaw.Dimension('x',-10.0,10.0,nx)
y = pyclaw.Dimension('y',-10.0,10.0,ny)
domain = pyclaw.Domain((x,y))

# Initialized state and Problem Parameters
state  = pyclaw.State(domain,3)
state.problem_data ={  'g':g , 'efix':True,'hr':True,'bcs':0,'cfix':1}
s_opts = {'f':0.1}

# Initial Data for 1d Dam Break
# Initial Data for 2d Radial Dam Break

h0fine = np.ones((1000,1000)) *H
xfine = pyclaw.Dimension('x',-10.0,10.0,1000)
yfine = pyclaw.Dimension('y',-10.0,10.0,1000)
rfine = np.sqrt(xfine.centers[:,None]**2 + yfine.centers[None,:]**2)
h0fine -= eta*np.sign(rfine - R)

inter = RectBivariateSpline(xfine.centers,yfine.centers,h0fine)

# Averaging the interpolant onto the coarse cells
h0 = np.zeros(state.grid.num_cells)
for i in xrange(x.num_cells):
    for j in xrange(y.num_cells):
        h0[i,j] = inter.integral(x.edges[i],x.edges[i+1],
                y.edges[j],y.edges[j+1])/np.prod(state.grid.delta)

r = np.sqrt(x.centers[:,None]**2 + y.centers[None,:]**2)
h0[r > 4] = H-eta

u0 = np.zeros(domain.grid.num_cells)
v0 = np.zeros(domain.grid.num_cells)

del xfine,yfine,h0fine,inter

# Figure out the time step
dx,dy = state.grid.delta
dt = min(dx,dy)/c/1.5
nt = int(T/dt)


state.q[0,:,:] = h0.copy()
state.q[1,:,:] = (u0*h0)
state.q[2,:,:] = (v0*h0)

cont = Controller(state,advance_sw)
cont.run(T)

from mpl_toolkits.mplot3d import Axes3D
rs = max(nx/100,1)
cs = max(ny/100,1)
xx,yy= state.grid.c_centers
fig = plt.figure(figsize=(12,4))
ax = fig.add_subplot(121,projection='3d')
ax.plot_surface(*cont.get_plot_args(),
    rstride=rs, cstride=cs,
    color='white',linewidth=.4,shade=True,alpha=1.0)

ax = fig.add_subplot(122,projection='3d')
ax.plot_surface(xx,yy,h0,
    rstride=rs, cstride=cs,
    color='white',linewidth=.4,shade=True,alpha=1.0)


# xx,yy,z = cont.get_plot_args()
# fig1 = plt.figure()
# plt.plot(xx[:,ny/2],z[:,ny/2])
# # plt.contourf(*cont.get_plot_args())
plt.show()


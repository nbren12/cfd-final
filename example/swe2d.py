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


if __name__=='__main__':

    g = .1
    H = 100
    nx = 100
    ny = 100
    ng = 2

    nframes = 10

    #   First Dimension is x, Second Dimension is y
    x = pyclaw.Dimension('x',-10.0,10.0,nx)
    y = pyclaw.Dimension('y',-5.0,5.0,ny)

    domain = pyclaw.Domain((x,y))
    state  = pyclaw.State(domain,3)
    state.problem_data ={  'g':.1 , 'efix':True,'hr':True,'bcs':0}
    s_opts = {'f':0.1}

    dx,dy = state.grid.delta
    dt = dx / 10
    T  = 1.0
    nt = int(T/dt)


    # Initial Data for 1D Riemann Problem
    # h0 = (2-np.sign(x.centers))[:,None] # Dam Break
    h0 = np.ones((nx,ny))*H #   + np.sign(x.centers)*H*.1
    # h0[nx/2-nx/3:nx/2+nx/3,ny/2-ny/3:ny/2+ny/3] = 3
    # h0[40:60,40:60] = 3
    # h0 += 2*np.exp(-reduce(lambda x,y:np.sqrt(x+y),map(lambda x: x**2,domain.grid.c_centers)))
    # h0 += 2*(reduce(lambda x,y:np.sqrt(x+y),map(lambda x: x**2,domain.grid.c_centers))< 4)
    u0 = np.zeros(domain.grid.num_cells)
    # u0 = +np.ones(domain.grid.num_cells)*np.sqrt(g*H)*.1
    u0[:,ny/2:ny/2+ny/6] = +np.sqrt(g*H)*.1
    u0[:,ny/2:] = np.sqrt(g*H)*2
    v0 = np.zeros(domain.grid.num_cells)


    xx,yy = state.grid.c_centers
    xe,ye = state.grid.c_edges
    state.q[0,:,:] = h0
    state.q[1,:,:] = (u0*h0)
    state.q[2,:,:] = (v0*h0)

    # cont = Controller(state,advance_sw,srcsolver=(advance_coriolis,s_opts))
    cont = Controller(state,advance_sw)

    for i in xrange(nt):
        cont.advance()

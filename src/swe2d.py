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

def periodic_bc2d(q,ng):
    nc,nx,ny = q.shape
    qbc = np.zeros(q.shape)

    qbc[:,ng:nx+ng,ng:ny+ng] = q

    return qbc
class Controller(object):
    """A Class to solve conservation laws"""

    def __init__(self,state,csolver,ssolver=None,num_ghost=2):
        """@todo: to be defined

        :state: @todo

        """
        self.state = state
        self.csolver = csolver

        if not None:
            self.ssolver = ssolver

        self.num_ghost = num_ghost

    def advance(self):
        """@todo: Docstring for advance
        :returns: @todo

        """
        dx,dy= state.grid.delta
        dt = min(dx/10,dy/10)
        advance_sw(self.state.q,dt,dx,dy,**state.problem_data)#   ,**self.state.problem_data)




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
state.problem_data ={  'g':.1 , 'efix':True,'hr':True}

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
cont = Controller(state,advance_sw)

for i in xrange(nt):
    # advance_sw(state.q,dt,dx,dy,**state.problem_data)
    cont.advance()

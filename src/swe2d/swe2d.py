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
import itertools
class Controller(object):
    """A Class to solve conservation laws"""

    def __init__(self,state,csolver,srcsolver=None):
        """@todo: to be defined

        :state: @todo

        """
        self.state = state
        self.csolver = csolver

        if srcsolver is not None:
            self.ssolver = srcsolver[0]
            self.s_opts = srcsolver[1]
        else:
            self.ssolver = None
            self.s_opts = None

        self.dt = min(state.grid.delta)/10


    def advance(self):
        """@todo: Docstring for advance
        :returns: @todo

        """
        delta = self.state.grid.delta
        self.csolver(self.state.q,self.dt,*delta,**self.state.problem_data)

        if self.ssolver is not None:
            self.ssolver(self.state.q,dt,**self.s_opts)
    def get_state(self):
        return self.state

    def run(self, T):
        """
        run model until time T

        :T: is the end time
        """
        nt = int(T/self.dt)

        for i in xrange(nt):
            print "Step %d of %d"%(i+1,nt)
            self.advance()
            self.time = ( i +1)*self.dt
            self.counter = i+1

    def get_plot_args(self,c=0):
        if self.state.grid.num_dim == 2:
            q = self.state.q[c,:,:]
        else:
            q = self.state.q[c,:]

        centers = self.state.grid.c_centers
        centers.append(q)
        return centers
    def get_plot_args_1d(self,i,c=0):
        x,y,z=  self.get_plot_args(c=c)
        pass

class SavedState(object):
    def __init__(self,name='swe',time=0.0):
        self.time = time
        self.name = name
        self.frame = -1

    def update_state(self, time,state):
        self.time  = time
        self.state = state
        self.frame +=1

    def init_dir(self):
        os.mkdir(self.name)


    def write_state(self):
        pass




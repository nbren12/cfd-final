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


    def advance(self):
        """@todo: Docstring for advance
        :returns: @todo

        """
        dx,dy= state.grid.delta
        dt = min(dx/10,dy/10)
        self.csolver(self.state.q,dt,dx,dy,**state.problem_data)

        if self.ssolver is not None:
            self.ssolver(self.state.q,dt,**self.s_opts)

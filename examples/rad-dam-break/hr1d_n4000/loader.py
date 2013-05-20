
import numpy as np
from clawpack import pyclaw
import os
from swe2d import Controller, ControllerSW2D,advance_sw
from scipy.interpolate import interp1d



# Initialize Grid
nx = 4000
ny = 2


#   First Dimension is x, Second Dimension is y
x = pyclaw.Dimension('x',0.0,15.0,nx)
y = pyclaw.Dimension('y',0.0,10.0,ny)
domain = pyclaw.Domain((x,y))

# Initialized state and Problem Parameters
state  = pyclaw.State(domain,3)
state.problem_data ={  'g':9.812000 , 'efix':1,'hr':1,'bcs':3,'cfix':0}

# Figure out the time step
dt = 0.000565
cont = ControllerSW2D(state,advance_sw,dt=dt,prefix=os.path.dirname(__file__))


def interp_maker(frame):
    cont.read_frame(frame)
    inter = interp1d(x.centers,cont.state.q[:,:,0])
    return inter

template = """
import numpy as np
from clawpack import pyclaw
import os
from swe2d import Controller, ControllerSW2D,advance_sw



# Initialize Grid
nx = %d
ny = %d


#   First Dimension is x, Second Dimension is y
x = pyclaw.Dimension('x',-10.0,10.0,nx)
y = pyclaw.Dimension('y',-10.0,10.0,ny)
domain = pyclaw.Domain((x,y))

# Initialized state and Problem Parameters
state  = pyclaw.State(domain,3)
state.problem_data ={  'g':%f , 'efix':%d,'hr':%d,'bcs':%d,'cfix':%d}

# Figure out the time step
dt = %f
cont = ControllerSW2D(state,advance_sw,dt=dt,prefix=os.path.dirname(__file__))
"""

def fill_template(nx=100,ny=100,g=9.812,bcs=0,efix=True,hr=True,cfix=True,dt=.02):
    return template%(nx,ny,g,efix,hr,bcs,cfix,dt)


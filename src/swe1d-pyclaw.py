"""
Shallow water Equations in One Dimension with no bottom topography are:

h_t + (hu)_x = 0
(hu)_t + (hu^2 + 1/2 g h^2 )_x = 0

hu^2 = q2^2/q1 + 1/2 g q2^2

Here, I used periodic boundary conditions. Might want to generalize the code later.
"""
import numpy as np
from clawpack import pyclaw
from clawpack.riemann.rp_shallow import rp_shallow_roe_1d

nx = 128
x = pyclaw.Dimension('x',-10.0,10.0,nx)
domain = pyclaw.Domain(x)
state = pyclaw.State(domain,2)

# Initial Data for 1D Riemann Problem
h0 = (2-np.sign(x.centers))
u0 = np.zeros(x.num_cells)

state.q[0,:]=h0
state.q[1,:]=u0*h0

state.problem_data['grav'] = 9.81
state.problem_data['efix'] = False

# Setup Solver
solver = pyclaw.ClawSolver1D()
solver.kernel_language = 'Python'
solver.rp = rp_shallow_roe_1d
solver.num_waves = 2
solver.bc_lower[0]=pyclaw.BC.wall
solver.bc_upper[0]=pyclaw.BC.wall


# Run solution
solution = pyclaw.Solution(state,domain)

controller = pyclaw.Controller()
controller.solution = solution
controller.solver   = solver
controller.tfinal   = 4
controller.run()
pyclaw.plot.html_plot()


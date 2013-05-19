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
import cPickle as pickle
import shutil
from .template import fill_template

def saveobject(obj, filename):
    with open(filename, 'wb') as output:
        pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)

def load_controller(folder):
    return pickle.load(folder)

class Controller(object):
    """A Class to solve conservation laws"""

    def __init__(self,state,csolver,srcsolver=None,dt=None,prefix='tmp_run'):
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

        if dt is not None:
            self.dt =dt
        else:
            self.dt = min(state.grid.delta)/10
        self.prefix = prefix
        self.folder = os.path.abspath(self.prefix)

    def save_controller(self):
        import __main__
        exists = os.path.exists(self.folder)
        if not exists: os.mkdir(self.folder)
        s = fill_template(
                nx=self.state.grid.num_cells[0],
                ny=self.state.grid.num_cells[1],
                dt=self.dt,
                **self.state.problem_data
                )
        with open(self.folder+"/loader.py",'w') as f:
            f.write(s)

        with open(self.folder+"/__init__.py",'a') as f:
            f.write('from .loader import cont')
        assert not exists
        shutil.copy2(os.path.realpath(__main__.__file__),self.folder)


        # saveobject(self,self.folder+"/controller.pickle")
    def save_frame(self, frame,time):
        """@todo: Docstring for save_frame

        :frame: @todo
        :time: @todo
        :returns: @todo

        """
        np.savez(self.folder+"/%.5d"%frame,t=time,q=self.state.q)

    def read_frame(self,frame):
        fpath = self.folder+"/%.5d.npz"%frame
        assert os.path.exists(fpath)
        npz =np.load(fpath)
        self.time = npz['t']
        self.state.q = npz['q']

    def advance(self):
        """@todo: Docstring for advance
        :returns: @todo

        """
        delta = self.state.grid.delta
        self.csolver(self.state.q,self.dt,*delta,**self.state.problem_data)

        if self.ssolver is not None:
            self.ssolver(self.state.q,self.dt,*delta,**self.s_opts)

    def get_state(self):
        return self.state

    def run(self, T,save=True):
        """
        run model until time T

        :T: is the end time
        a = qb +r
        """
        fac = int(0.1/self.dt)
        self.dt= .1/fac

        nt = int(T/self.dt)

        if save:
            self.save_controller()
            self.save_frame(0,0)
            self.frame=0
            dframe = fac

        for i in xrange(nt):
            self.time = ( i +1)*self.dt
            print "Step %d of %d"%(i+1,nt)
            self.advance()
            if save and (i+1)%dframe == 0:
                print "Saving Frame at Time = %f"%self.time
                self.frame +=1
                self.save_frame(self.frame,self.time)

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

class ControllerSW2D(Controller):
    def surf_plot(self,ax=None):
        nx,ny = self.state.grid.num_cells
        rs = max(nx/100,1)
        cs = max(ny/100,1)
        ax.plot_surface(*self.get_plot_args(),
            rstride=rs, cstride=cs,
            color='white',linewidth=.4,shade=True,alpha=1.0)
    def cont_plot(self,ax=None):
        xx,yy,h = self.get_plot_args()
        ax.contour(xx,yy,h,20)




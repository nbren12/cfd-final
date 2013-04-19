import numpy as np
"""
Shallow water Equations in One Dimension with no bottom topography are:

h_t + (hu)_x = 0
(hu)_t + (hu^2 + 1/2 g h^2 )_x = 0


Here, I used periodic boundary conditions. Might want to generalize the code later.
"""

class OneDGrid:
    def __init__(self,x=None,m=None,bounds=[-1.0,1.0]):
        if m is not None:
            m = float(m)
            L = bounds[1]-bounds[0]
            h = L/m
            x = h*(np.arange(m) + .5) + bounds[0]
            xl = np.arange(m)*h + bounds[0]
        elif x is not None:
            m = float(x.shape[0])
            h = L/m

        self.bounds = bounds
        self.L = L
        self.m = int(m)
        self.h = h
        self.x = x
        self.xl = xl

class SWSolution:
    def __init__(self,grid=None):
        """docstring for __init__"""
        pass

def ReimannSolve(ql,qr,A):
    pass

def CalcFlux(ul,ur,hl,hr):
    return (1,1)


if __name__=='__main__':
    gg = OneDGrid(m=128)
    ul = np.zeros(gg.m)
    hl = np.zeros(gg.m)
    u0 = np.zeros(gg.m)
    h0 = (np.sign(gg.x) -1) / (-2) + 1
    u = np.zeros(gg.m)
    h = np.zeros(gg.m)

    # Setup Time Stepping
    tau = gg.h/2
    T   = 1.0
    nt  = T/tau + 1
    t   = np.arange(nt)*tau

    for tt in np.nditer(t[1:]):

        # Calculate Fluxes
        for i in range(gg.m):
            ul[i],hl[i] = CalcFlux(u[i-1],u[i],h[i-1],h[i])

from pylab import *
import glob
import os
from scipy.interpolate import interp1d
from ipdb import set_trace as st
frame = 6
suffix = "o1"
ion()
def calc_error(cont1,cont2):

    cont1.read_frame(frame)
    cont2.read_frame(frame)

    q1 = cont1.state.q
    h1 = q1[0,:,0]
    x1 = cont1.state.grid.x.centers

    #    periodixc
    x11 = cont1.state.grid.x.centers_with_ghost(1)
    h11 = np.zeros(h1.shape[0] + 2)
    h11[1:-1] = h1
    h11[0] = h1[-1]
    h11[-1] = h1[0]


    q2 = cont2.state.q
    h2 = q2[0,:,0]
    x2 = cont2.state.grid.x.centers

    inter  = interp1d(x11,h11,'nearest')

    h1_up = inter(x2)
    err2 = ((h1_up-h2)**2).sum()*cont2.state.grid.delta[0]
    return sqrt(err2)

p = []
ll =  [2**n for n in range(5,12)]
figure()
for i in xrange(len(ll)-1):
    m = __import__("geo_break1d_n%d_%s"%( ll[i],suffix ))
    m1 = __import__("geo_break1d_n%d_%s"%(ll[i+1],suffix))

    err = calc_error(m.cont,m1.cont)
    print err
    p.append(err)


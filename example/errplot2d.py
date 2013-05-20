from pylab import *
import glob
import os
from ipdb import set_trace as st

from hr1d_n4000 import interp_maker
from hr1d_n4000 import cont as hrdat

datadirs = glob.glob('n*_2dhr')
datadirs = filter(
        lambda x : os.path.isdir(x),
        datadirs)

frame= 7
hrcont = interp_maker(frame)
hrdat.read_frame(frame)

def calc_error(cont):
    q = cont.state.q
    h = q[0,:,:]
    u = q[1,:,:]/h
    v = q[2,:,:]/h

    ur = sqrt(u**2 + v**2)


    x,y = cont.state.grid.c_centers

    r = sqrt(x**2 + y**2)
    exact = hrcont(r)
    he = exact[0,:,:]


    err2 = ((he-h)**2).sum()/(he**2).sum()

    return sqrt(err2),r,q


for dd in datadirs:
    print dd
    data = __import__(dd)
    data.cont.read_frame(frame)
    assert(isnan(data.cont.state.q).sum()==0)



    # err= calc_error(data.cont)
    cont = data.cont
    err , r,q = calc_error(cont)
    print err

figure()
plot(r.ravel(),q[0,:,:].ravel(),'o')
plot(hrdat.state.grid.x.centers,hrdat.state.q[0,:,0])
show()


import numpy as np


def advance_metric_terms(q,dt,dx,dy,r=None):
    qm = q.copy()
    qm[0,:,:] = q[0,:,:] - (dt/2.0)*q[0,:,:]*q[1,:,:]/r
    qm[1,:,:] = q[1,:,:] - (dt/2.0)*q[0,:,:]*q[1,:,:]**2/r

    q[0,:,:] = q[0,:,:] - (dt/2.0)*qm[0,:,:]*qm[1,:,:]/r
    q[1,:,:] = q[1,:,:] - (dt/2.0)*qm[0,:,:]*qm[1,:,:]**2/r

    return q

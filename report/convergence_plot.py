from pylab import *


s = genfromtxt('convergence2d.csv',delimiter=',',skip_header=1)


n = s[:,0]
h = 20.0 / n
hr = s[:,1]
o1 = s[:,2]
hr_nocfix = s[:,3]
hr_noefix = s[:,4]
fig = figure()

ax =fig.add_subplot(111)

ax.loglog(h,hr,'ko-')
ax.loglog(h,o1,'kx-')
ax.loglog(h,hr_nocfix,'bD-')
ax.loglog(h,hr_noefix,'r^-')

ps = []
for ll in [hr,o1,hr_nocfix,hr_noefix]:
    p = polyfit(log(h),log(ll),1)
    ps.append(p[0])

ss = ['HR w/ cfix','Godunov w/ cfix','HR w/o cfix','HR w/o efix']
s = [ss[i] + " p=%.2f"%ps[i] for i in range(4)]

ax.set_xlim([h[-1],h[0]])
ax.set_ylim([5e-3,1e-1])
ax.legend(s,loc=4)
ax.grid(b=True,which='major',linestyle='-')
ax.grid(b=True,which='minor',linestyle='--')

ax.set_xlabel('Step Size')
ax.set_ylabel('Relative L2 Error')
fig.savefig('convergence_raddam_2d.eps',bbox_inches=None)

show()

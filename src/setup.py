#/usr/bin/env python

from numpy.distutils.core import setup, Extension
from numpy.distutils.misc_util import Configuration

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration(None,parent_package,top_path)


    config.add_extension('solver2d',
                    ['./hr_solver2d.f90','./rp_sw2d_roe.f90','./bc2d.f90'],
                    )


    return config



setup(configuration=configuration)

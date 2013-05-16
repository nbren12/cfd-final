#/usr/bin/env python

from numpy.distutils.core import setup, Extension
from numpy.distutils.misc_util import Configuration

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('swe2d',parent_package,top_path)


    config.add_extension('solver',
                    ['./hr_solver2d.f90','./rp_roe.f90','./bc2d.f90'],
                    )


    return config


if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(configuration=configuration)

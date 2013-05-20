#/usr/bin/env python

from numpy.distutils.core import setup, Extension
from numpy.distutils.misc_util import Configuration

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('',parent_package,top_path)


    config.add_extension('adv_weno',
                    ['./adv_weno.f90','./weno.f90'],
                    )


    return config


if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(configuration=configuration)

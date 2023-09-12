from setuptools import setup

from hadronic_module import __version__

setup(
    name='hadronic_module',
    version=__version__,

    url='https://github.com/micro-uhecr/prototype_him_crpropa',
    author='Leonel Morejon',
    author_email='leonel.morejon@uni-wuppertal.de',

    py_modules=['hadronic_module'],
)
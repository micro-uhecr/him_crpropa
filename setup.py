from setuptools import setup

from him_crpropa.hadronic_module import __version__

setup(
    name='him_crpropa',
    version=__version__,

    url='https://github.com/micro-uhecr/prototype_him_crpropa',
    author='Leonel Morejon',
    author_email='leonel.morejon@uni-wuppertal.de',
    
    packages=['him_crpropa'],
    install_requires=['numpy','scipy', 'chromo>=0.5.0'],
    py_modules=['him_crpropa'],
)
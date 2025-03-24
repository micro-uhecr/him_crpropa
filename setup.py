from setuptools import setup

setup(
    name='him_crpropa',
    version="0.1.0",

    url='https://github.com/micro-uhecr/prototype_him_crpropa',
    author='Leonel Morejon',
    author_email='leonel.morejon@uni-wuppertal.de',
    
    install_requires=['numpy','scipy', 'chromo>=0.5.0'],
    packages=['him_crpropa'],
    py_modules=['him_crpropa'],
)
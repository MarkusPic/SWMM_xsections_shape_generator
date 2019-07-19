__author__ = "Markus Pichler"
__credits__ = ["Markus Pichler"]
__maintainer__ = "Markus Pichler"
__email__ = "markus.pichler@tugraz.at"
__version__ = "0.2"
__license__ = "MIT"

from setuptools import setup

setup(
    name='SWMM_xsections_shape_generator',
    version='0.2',
    packages=['shape_generator'],
    url='https://github.com/tugraz-sww/SWMM_xsections_shape_generator',
    license='MIT',
    author='Markus Pichler',
    author_email='markus.pichler@tugraz.at',
    description='US-EPA SWMM Cross-Section curve shape generator',  # TODO
    install_requires=['numpy', 'pandas', 'matplotlib', 'sympy'],
)

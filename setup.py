__author__ = "Markus Pichler"
__credits__ = ["Markus Pichler"]
__maintainer__ = "Markus Pichler"
__email__ = "markus.pichler@tugraz.at"
__version__ = "0.1"
__license__ = "MIT"

from setuptools import setup

setup(
    name='SWMM_xsections_shape_generator',
    version='0.1',
    packages=['shape_generator'],
    url='https://github.com/MarkusPic/SWMM_xsections_shape_generator',
    license='MIT',
    author='Markus Pichler',
    author_email='markus.pichler@tugraz.at',
    description='Diverse tools to export and analyse the >10a rain series from the ehyd.gv.at platform',
    install_requires=['numpy', 'pandas', 'matplotlib', 'sympy'],
)

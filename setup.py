__author__ = "Markus Pichler"
__credits__ = ["Markus Pichler"]
__maintainer__ = "Markus Pichler"
__email__ = "markus.pichler@tugraz.at"
__version__ = "0.2"
__license__ = "MIT"

from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()

with open('requirements.txt') as f:
    requirements = f.read().splitlines()

setup(
    name='SWMM_xsections_shape_generator',
    version='0.3.2',
    packages=['shape_generator'],
    url='https://github.com/MarkusPic/SWMM_xsections_shape_generator',
    license='MIT',
    author='Markus Pichler',
    author_email='markus.pichler@tugraz.at',
    description='US-EPA SWMM Cross-Section curve shape generator',  # TODO
    install_requires=requirements,
    long_description=long_description,
    long_description_content_type="text/markdown",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        'Topic :: Scientific/Engineering',
        'Topic :: Software Development :: Libraries',
        'Topic :: Software Development :: Libraries :: Application Frameworks',
        # 'Framework :: SWMM',
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Natural Language :: English'
    ],
    keywords='swmm environment civil engineering api',
    python_requires='>= 3.5',
)

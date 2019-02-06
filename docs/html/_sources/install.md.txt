
# Install

The script is written in Python3.

## Windows

I recommend to use [Anaconda](https://www.anaconda.com/download/) to install python on Windows and the Anaconda-Prompt for the commandline tool.

Alternatively, you can install the original python from the [website](https://www.python.org/downloads/).
To use the syntax explained in the usage section below, 
you have to add the path to your python binary to the environment variables. 
This is an option in the installation window as seen below:

- [x] Add Python 3.7 to PATH

![python_install](https://github.com/MarkusPic/ehyd_tools/raw/master/example/python_install.png)

## Linux/Unix

Python is pre-installed on most operating systems.

## python Packages

Packages required for this program will be installed with pip during the installation process and can be seen in the 'requirements.txt' file.

## Fresh install


```
pip install https://codeload.github.com/MarkusPic/SWMM_xsections_shape_generator/zip/master
```

To install the package only for the local user account, add ```--user``` to the install command.

## Update package

To update the package, add ```--upgrade``` to the install command.

```
pip install https://codeload.github.com/MarkusPic/SWMM_xsections_shape_generator/zip/master --upgrade
```

from sympy import Expr, sqrt, Symbol, solve, Float
from pandas import read_csv, DataFrame
import pandas
import sympy
from io import StringIO
from os import path, listdir, remove
from math import radians, cos, tan

x = Symbol('x', real=True, positive=True)
d = Symbol('d', real=True)


def csv(txt, comment=None):
    """Read the string in txt as csv file

    Read the string in txt as csv file and return the content as DataFrame.

    Args:
        txt (str): content of csv
        comment (str): comment sign

    Returns:
        pandas.DataFrame: csv table as pandas DataFrame
    """
    df = read_csv(StringIO(txt), index_col=0, skipinitialspace=True, skip_blank_lines=True, comment=comment)
    df = df[df.index.notnull()].copy()
    df.index = df.index.astype(str)
    return df


def deg2slope(deg):
    """

    convert degrees to a slope (:math:`\\Delta y / \\Delta x`)

    Args:
        deg (float): angle in degree

    Returns:
        float: slope
    """
    return round(tan(radians(deg)), 5)


def channel_end(r, alpha):
    """
    get vertical end of the channel

    Args:
        r (float):
        alpha (float): in degree

    Returns:
        float: height of the channel when the circle reaches a certain angle
    """
    return r * (1 - cos(radians(alpha)))


def circle(r, x_m=0, y_m=0, clockwise=False):
    """
    get funktion of a circle with a given mid point

    Args:
        r (float):  radius
        x_m (float): x axis value of the mid point
        y_m (float): y axis value of the mid point
        clockwise (bool): whether the circle is clockwise or anticlockwise

    Returns:
        sympy.core.expr.Expr: function of the circle
    """
    return sqrt(Float(float(r)) ** 2 - (x - Float(float(x_m))) ** 2) * (-1 if clockwise else 1) + Float(float(y_m))


def linear(slope, p0):
    """

    Args:
        slope (float):
        p0 (float):

    Returns:
        sympy.core.expr.Expr: linear function
    """
    x0, y0 = p0
    if slope == 0:
        return x0
    fi = (x - d) / slope
    di = solve(fi.subs(x, x0) - y0, d)[0]
    fi = fi.subs(d, di)
    return fi


def combine_input_files(shape_path, delete_original=False):
    """combine all generated shape text files to a single inp-like text file

    When running the :func:`shape_generator.shape_generator.Profile.input_file` function, a .txt file will be created.
    Those txt files will be combines to a single file with this function.
    This makes it easier to import all shapes to the .inp file.

    Args:
        shape_path (str): path where the shapes are stored
        delete_original (bool): whether to delete the original single files
    """
    with open(path.join(shape_path, 'all_shapes.txt'), 'w') as outfile:
        for fname in listdir(shape_path):
            if not fname.endswith('_shape.txt'):
                continue
            in_fn = path.join(shape_path, fname)
            with open(in_fn) as infile:
                outfile.write('\n')
                outfile.write(infile.read())
            if delete_original:
                remove(in_fn)
    print('Files are combined and originals {}deleted.'.format('' if delete_original else 'NOT '))

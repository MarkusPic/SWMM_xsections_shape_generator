from sympy import Expr, sqrt, Symbol, solve, Float
from pandas import read_csv, DataFrame
from io import StringIO
from os import path, listdir, remove
from math import radians, cos, tan

x = Symbol('x', real=True, positive=True)
d = Symbol('d', real=True)


def csv(txt, comment=None):
    """
    Read the string in txt as csv file and return the content as DataFrame.

    :param txt: content of csv
    :type txt: str

    :param comment: comment sign
    :type comment: str

    :return: csv table as pandas DataFrame
    :rtype: DataFrame
    """
    df = read_csv(StringIO(txt), index_col=0, skipinitialspace=True, skip_blank_lines=True, comment=comment)
    df = df[df.index.notnull()].copy()
    df.index = df.index.astype(str)
    return df


def deg2slope(deg):
    """
    convert degrees to a slope (delta y / delta x)

    :math:`\\delta`

    :param deg: angle in degree
    :type deg: float

    :return: slope
    :rtype: float
    """
    return round(tan(radians(deg)), 5)


def channel_end(r, alpha):
    """
    get vertical end of the channel

    :param float r:

    :param float alpha: in grad

    :return:
    :rtype: float
    """
    return r * (1 - cos(radians(alpha)))


def circle(r, x_m=0, y_m=0, uzs=False):
    """

    :param float r:
    :param float x_m:
    :param float y_m:
    :param bool uzs:
    :return: circle function
    :rtype: Expr
    """
    return sqrt(Float(float(r)) ** 2 - (x - Float(float(x_m))) ** 2) * (-1 if uzs else 1) + Float(float(y_m))


def linear(slope, p0):
    """

    :param float slope:

    :param tuple p0:

    :return: linear function
    :rtype: Expr
    """
    x0, y0 = p0
    if slope == 0:
        return x0
    fi = (x - d) / slope
    di = solve(fi.subs(x, x0) - y0, d)[0]
    fi = fi.subs(d, di)
    return fi


def combine_input_files(shape_path, delete_original=False):
    """

    :param shape_path: path where the shapes are stored

    :param delete_original:

    :return:
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

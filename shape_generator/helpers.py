from io import StringIO
from os import path, listdir, remove
from math import radians, tan, cos, pi, atan, sin

import pandas as pd
import sympy as sy

# these variables are used to solve symbolic mathematical equations
# x is the control variable over the height ... max(x) = H_cross_section
x = sy.Symbol('x', real=True, positive=True)
d = sy.Symbol('d', real=True)  # interception of the linear function


def csv(txt, comment=None):
    """
    Read the string in txt as csv file and return the content as DataFrame.

    Args:
        txt (str): content of csv
        comment (str): comment sign

    Returns:
        pandas.DataFrame: csv table as pandas DataFrame
    """
    df = pd.read_csv(StringIO(txt), index_col=0, skipinitialspace=True, skip_blank_lines=True, comment=comment)
    df = df[df.index.notnull()].copy()
    df.index = df.index.astype(str)
    return df


def deg2slope(degree):
    """
    convert degrees to a slope (:math:`\\Delta x / \\Delta y`)

    Args:
        degree (float): angle in degree

    Returns:
        float: slope

    .. figure:: images/slope.gif
        :align: center
        :alt: slope
        :figclass: align-center

        Slope
    """
    return round(tan(radians(degree)), 5)


def channel_end(r, end_degree):
    """
    get vertical end of the channel based on the radius of the channel and an end angle

    Args:
        r (float): radius of the channel
        end_degree (float): end angle in degree (°)

    Returns:
        float: height of the channel when the circle reaches a certain angle

    .. figure:: images/channel_end.gif
        :align: center
        :alt: channel end
        :figclass: align-center

        Channel end
    """
    return r * (1 - cos(radians(end_degree)))


# def linear(slope, p0):
#     """
#     get function/expression of a straight line with a given point which it intercepts
#
#     Args:
#         slope (float): slope
#         p0 (set[float, float]): point as a set of a x and a y coordinate
#
#     Returns:
#         sympy.core.expr.Expr: linear function
#
#     .. figure:: images/gerade.gif
#         :align: center
#         :alt: straight line
#         :figclass: align-center
#
#         Straight line
#     """
#     x0, y0 = p0
#     if slope == 0:
#         return x0
#     fi = (x - d) / slope
#     di = sy.solve(fi.subs(x, x0) - y0, d)[0]
#     fi = fi.subs(d, di)
#     return fi


def sqrt(i):
    return i ** (1 / 2)


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


####################################################################################################################
# testing new functions
####################################################################################################################
# def interp(var_x, x0, x1, f0, f1):
#     """
#     create sympy function for linear interpolation between two points
#
#     Args:
#         var_x (sympy.Symbol): x
#         x0 (float):
#         x1 (float):
#         f0 (float):
#         f1 (float):
#
#     Returns:
#         sympy.Expr: linear interpolation between two points
#     """
#     return f0 + (f1 - f0) / (x1 - x0) * (var_x - x0)
#
#
# def solve_equation(f, xi):
#     """
#     solve a given equation
#
#     Args:
#         f (sympy.Expr):
#         xi (float):
#
#     Returns:
#         float: result
#     """
#     if isinstance(f, sy.Expr):
#         return float(f.subs(x, sy.Float(round(xi, 3))))
#     elif isinstance(f, CustomExpr):
#         return


class CustomExpr:
    def __init__(self):
        pass

    def __repr__(self):
        return 'Custom Function'

    def expr(self):
        pass

    def solve(self, i):
        pass

    def length(self, i0, i1):
        pass

    def area(self, i0, i1):
        pass


class Slope(CustomExpr):
    def __init__(self, slope, unit=None):
        if unit is None:
            self.slope = slope
        elif unit == '°':
            self.slope = deg2slope(slope)
        elif unit == '%':
            self.slope = slope / 100

        self.x0 = None
        self.y0 = None

        self.x1 = None
        self.y1 = None

        CustomExpr.__init__(self)

    def __repr__(self):
        return f'Slope Function (k={self.slope:0.2f}, zero=[{self.x0:0.2f}, {self.y0:0.2f}])'

    def set_start_point(self, x0, y0):
        self.x0 = x0
        self.y0 = y0

    def set_end_point(self, x1, y1):
        self.x1 = x1
        self.y1 = y1

    def expr(self):
        return self.y0 + (x - self.x0) / self.slope

    def solve(self, i):
        return self.y0 + (i - self.x0) / self.slope

    @staticmethod
    def from_points(start, end):
        x0, f0 = start
        x1, f1 = end
        slope = (f1 - f0) / (x1 - x0)
        new_slope = Slope(slope)
        new_slope.set_start_point(x0, f0)
        new_slope.set_start_point(x1, f1)
        return new_slope

    def end_point(self):
        pass

    def length(self, i0, i1):
        return sqrt((self.solve(i0) - self.solve(i1)) ** 2 + (i0 - i1) ** 2)

    def area(self, i0, i1):
        return (self.solve(i0) + self.solve(i1)) / 2 * abs(i0 - i1)


class Vertical(CustomExpr):
    def __init__(self, y):
        self.y = y
        CustomExpr.__init__(self)

    def __repr__(self):
        return f'Vertical Function (y={self.y:0.2f})'

    def expr(self):
        return self.y + x * 0

    def solve(self, i):
        return self.y + i * 0

    def length(self, i0, i1):
        return i1 - i0

    def area(self, i0, i1):
        return self.length(i0, i1) * self.y


def d_alpha(xm, ym, h, b):
    if (b - ym) == 0:
        a = pi / 2
        if (h - xm) < 0:
            a *= -1
    else:
        a = atan((h - xm) / (b - ym))
    return a


class Circle(CustomExpr):
    def __init__(self, r, x_m=0, y_m=0, clockwise=False):
        self.r = float(r)
        self.x_m = float(x_m)
        self.y_m = float(y_m)
        self.clockwise = clockwise

        CustomExpr.__init__(self)

    def __repr__(self):
        return f'Circle Function (radius={self.r:0.2f}, mid=[{self.x_m:0.2f}, {self.y_m:0.2f}])'

    def expr(self):
        """
        get function/expression of a circle with a given mid point

        Args:
            r (float): radius
            x_m (float): x axis value of the mid point
            y_m (float): y axis value of the mid point
            clockwise (bool): whether the circle is clockwise or anticlockwise

        Returns:
            sympy.core.expr.Expr: function of the circle

        .. figure:: images/kreis.gif
            :align: center
            :alt: circle
            :figclass: align-center

            Circle
        """
        return sy.sqrt(sy.Float(self.r) ** 2 - (x - sy.Float(self.x_m)) ** 2) * (-1 if self.clockwise else 1) + \
               sy.Float(self.y_m)

    def solve(self, i):
        return sqrt(self.r ** 2 - (i - self.x_m) ** 2) * (-1 if self.clockwise else 1) + self.y_m

    def length(self, i0, i1):
        alpha = abs(d_alpha(self.x_m, self.y_m, i0, self.solve(i0)) - d_alpha(self.x_m, self.y_m, i1, self.solve(i1)))
        return alpha * self.r

    def area(self, i0, i1):
        alpha = abs(d_alpha(self.x_m, self.y_m, i0, self.solve(i0)) - d_alpha(self.x_m, self.y_m, i1, self.solve(i1)))
        return self.r ** 2 / 2 * (alpha - sin(alpha)) + (self.solve(i0) + self.solve(i1)) / 2 * self.length(i0, i1)

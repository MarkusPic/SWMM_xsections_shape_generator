import io
import math
import os
import warnings

import numpy as np
import pandas as pd

accuracy = 10


def to_num(x):
    if x == '':
        return None
    elif x.replace('-', '').isdecimal():
        return int(x)
    elif ('.' in x) and (x.lower().replace('.', '').replace('-', '').replace('e', '').isdecimal()):
        return float(x)
    else:
        return x


def csv(txt, comment=None):
    """
    Read the string in txt as csv file and return the content as DataFrame.

    Args:
        txt (str): content of csv
        comment (str): comment sign

    Returns:
        dict: profile label and values
    """
    df = pd.read_csv(io.StringIO(txt), index_col=0, skipinitialspace=True, skip_blank_lines=True, comment=comment)
    df = df[df.index.notnull()].copy()
    df.index = df.index.astype(str)
    return df


def to_xs_dict(txt, comment=None):
    """
    Read the string in txt as csv file and return the content as DataFrame.

    Args:
        txt (str): content of csv
        comment (str): comment sign

    Returns:
        dict: profile label and values
    """
    di = dict()
    names = []
    for line in txt.split('\n'):
        if line == '':
            continue
        elif isinstance(comment, str) and line.startswith(comment):
            continue
        elif not names:
            names = [n.strip() for n in line.split(',')[1:]]
            di['_names'] = names
        else:
            name, *values = [n.strip() for n in line.split(',')]
            # di[name] = {k: to_num(v) for k, v in zip(names, values)}
            di[name] = [to_num(v) for v in values]
    return di


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
    return math.tan(math.radians(degree))


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
    return r * (1 - math.cos(math.radians(end_degree)))


def combine_input_files(shape_path, delete_original=False):
    """combine all generated shape text files to a single inp-like text file

    When running the :func:`shape_generator.shape_generator.Profile.input_file` function, a .txt file will be created.
    Those txt files will be combines to a single file with this function.
    This makes it easier to import all shapes to the .inp file.

    Args:
        shape_path (str): path where the shapes are stored
        delete_original (bool): whether to delete the original single files
    """
    with open(os.path.join(shape_path, 'all_shapes.txt'), 'w') as outfile:
        for fname in os.listdir(shape_path):
            if not fname.endswith('_shape.txt'):
                continue
            in_fn = os.path.join(shape_path, fname)
            with open(in_fn) as infile:
                outfile.write(infile.read())
                outfile.write('\n\n')
            if delete_original:
                os.remove(in_fn)
    print(f'Files are combined and originals{" " if delete_original else " NOT "}deleted.')


def get_intersection_point(expr1, expr2, x_from, x_to):
    from scipy.optimize import minimize_scalar
    x_i = minimize_scalar(lambda j: abs(expr1.solve_y(j) - expr2.solve_y(j)), bounds=(x_from, x_to), method='bounded', options=dict(xatol=.01)).x
    # y_i = float(expr1.solve_y(x_i))

    # import sympy as sy
    # # these variables are used to solve symbolic mathematical equations
    # # x is the control variable over the height ... max(x) = H_cross_section
    # x = sy.Symbol('x', real=True, positive=True)
    #
    # expr_eq = expr1.expr(x) - expr2.expr(x)
    #
    # res = sy.solve(expr_eq, x)
    # if len(res) == 0:
    #
    #     x_i = minimize_scalar(lambda j: float(expr_eq.subs(x, j)),
    #                           bounds=(x_from, x_to), method='bounded').x
    #
    # elif len(res) == 1:
    #     x_i = float(res[0])
    # else:
    #     # multiple results
    #     # TODO: how to handle it
    #     x_i = float(res[0])

    y_i = float(expr1.solve_y(x_i))
    return x_i, y_i


####################################################################################################################
class CustomExpr:
    def __init__(self):
        self.x0 = None
        self.y0 = None

        self.x1 = None
        self.y1 = None

    def __repr__(self):
        return 'Custom Function'

    def expr(self, x):
        pass

    def solve_y(self, i):
        pass

    def solve_x(self, y):
        pass

    def length(self, i0, i1):
        pass

    def area(self, i0, i1):
        pass

    def _fix_point(self, point):
        x, y = point
        if x is None:
            x = float(self.solve_x(y))

        if y is None:
            y = float(self.solve_y(x))
        return x, y

    def set_start_point(self, point):
        """set start point"""
        x0, y0 = self._fix_point(point)
        self.x0 = x0
        self.y0 = y0

    def set_end_point(self, point):
        """set end point"""
        x1, y1 = self._fix_point(point)
        self.x1 = x1
        self.y1 = y1

        if self.get_start_point() == self.get_end_point():
            warnings.warn('unused part of the shape detected. Ignoring this part.')

    def get_start_point(self):
        """get the start point"""
        return self.x0, self.y0

    def get_end_point(self):
        """get the end point"""
        return self.x1, self.y1


####################################################################################################################
class Linear(CustomExpr):
    def __init__(self):
        CustomExpr.__init__(self)

    def set_points(self, start, end):
        self.set_start_point(start)
        self.set_end_point(end)

    @classmethod
    def from_points(cls, start, end):
        """
        set the slope by giving the start and end point

        Args:
            start ():
            end ():

        Returns:

        """
        x0, y0 = start
        x1, y1 = end

        if abs(y0 - y1) < 1.0e-6:
            new_obj = Vertical(y0)
        elif abs(x0 - x1) < 1.0e-6:
            new_obj = Horizontal(x0)
        else:
            slope = (x1 - x0) / (y1 - y0)
            new_obj = Slope(slope)

        new_obj.set_points(start, end)
        return new_obj


class Slope(Linear):
    """
    get function/expression of a straight line with a given point which it intercepts

    Args:
        slope (float): slope
        p0 (set[float, float]): point as a set of a x and a y coordinate

    Returns:
        sympy.core.expr.Expr: linear function

    .. figure:: images/gerade.gif
        :align: center
        :alt: straight line
        :figclass: align-center

        Straight line
    """

    def __init__(self, slope, unit=None):
        if unit is None or unit == '':
            self.slope = slope
        elif unit == '°':
            self.slope = deg2slope(slope)
        elif unit == '%':
            self.slope = slope / 100
        else:
            raise NotImplementedError('Unknown Unit for slope function')

        Linear.__init__(self)

    def __repr__(self):
        return f'Slope Function (k={self.slope:0.2f}, zero=[{self.x0:0.2f}, {self.y0:0.2f}])'

    def expr(self, x):
        """get sympy expression"""
        return self.y0 + (x - self.x0) / self.slope

    def solve_y(self, x):
        """get y value"""
        return self.y0 + (x - self.x0) / self.slope

    def solve_x(self, y):
        return self.x0 + (y - self.y0) * self.slope

    def length(self, i0, i1):
        """get shape length between two values"""
        return math.sqrt((self.solve_y(i0) - self.solve_y(i1)) ** 2 + (i0 - i1) ** 2)

    def area(self, x0, x1):
        """get shape area between two values"""
        return (self.solve_y(x0) + self.solve_y(x1)) / 2 * np.abs(x0 - x1)


####################################################################################################################
class Vertical(Linear):
    """
    function of a vertical line

    for a shape curve this means a constant width
    """
    def __init__(self, y=None):
        """

        Args:
            y (float): y value of the vertical line
        """
        Linear.__init__(self)
        if y is not None:
            self.set_y(y)

    def __repr__(self):
        return f'Vertical Function (y={self.y:0.2f})'

    def expr(self, x):
        return self.y + x * 0

    def solve_y(self, x):
        if isinstance(x, (list, tuple, np.ndarray)):
            return [self.y]*len(x)
        else:
            return self.y

    def solve_x(self, y):
        pass

    def length(self, x0, x1):
        return x1 - x0

    def area(self, x0, x1):
        return self.length(x0, x1) * self.y

    @property
    def y(self):
        return self.y0

    def set_y(self, y):
        self.y0 = y
        self.y1 = y


####################################################################################################################
class Horizontal(Linear):
    """
    function of a horizontal line
    """
    def __init__(self, x=None):
        Linear.__init__(self)
        if x is not None:
            self.set_x(x)

    def __repr__(self):
        return 'Horizontal Function'

    def expr(self, x):
        return self.y1

    def solve_y(self, i):
        pass

    def solve_x(self, y):
        if isinstance(y, (list, tuple, np.ndarray)):
            return [self.x]*len(y)
        else:
            return self.x

    def length(self, i0, i1):
        return np.abs(self.y1 - self.y0)

    def area(self, i0, i1):
        return 0

    @property
    def x(self):
        return self.x0

    def set_x(self, x):
        self.x0 = x
        self.x1 = x


####################################################################################################################
class Circle(CustomExpr):
    """
    function of a circle

    .. figure:: images/kreis.gif
        :align: center
        :alt: circle
        :figclass: align-center

        Circle
    """
    def __init__(self, r, x_m=0, y_m=0, clockwise=False):
        """

        Args:
            r (float): radius
            x_m (float): x axis value of the mid point
            y_m (float): y axis value of the mid point
            clockwise (bool): whether the circle is clockwise or anticlockwise

        """

        self.r = float(r)
        self.x_m = float(x_m)
        self.y_m = float(y_m)
        self.clockwise = clockwise

        CustomExpr.__init__(self)

    def __repr__(self):
        return f'Circle Function (radius={self.r:0.2f}, mid=[{self.x_m:0.2f}, {self.y_m:0.2f}])'

    def expr(self, x):
        """
        get function/expression of a circle with a given mid point


        Returns:
            sympy.core.expr.Expr: function of the circle
        """
        import sympy as sy
        return sy.sqrt(sy.Float(self.r) ** 2 - (x - sy.Float(self.x_m)) ** 2) * (-1 if self.clockwise else 1) + \
               sy.Float(self.y_m)

    def _alpha(self, i):
        """
        angle in the circle of a point to the horizontal

        Args:
            i: variable

        Returns:
            float: angle in rad
        """
        if isinstance(i, np.ndarray):
            return np.arctan((i - self.x_m) / (self.solve_y(i) - self.y_m))

        else:
            if (self.solve_y(i) - self.y_m) == 0:
                a = math.pi / 2
                if (i - self.x_m) < 0:
                    a *= -1
            else:
                a = np.arctan((i - self.x_m) / (self.solve_y(i) - self.y_m))
            return a

    def _d_alpha(self, i0, i1):
        """
        difference of the angle in the circle of two points

        Args:
            i0: start variable
            i1: end variable

        Returns:
            float: difference of the angle in rad
        """
        return np.abs(self._alpha(i0) - self._alpha(i1))

    def solve_y(self, x):
        return np.sqrt(self.r ** 2 - (x - self.x_m) ** 2) * (-1 if self.clockwise else 1) + self.y_m

    def solve_x(self, y):
        return np.sqrt(self.r ** 2 - (y - self.y_m) ** 2) * (-1 if self.clockwise else 1) + self.x_m

    def length(self, i0, i1):
        return self._d_alpha(i0, i1) * self.r

    def area(self, i0, i1):
        alpha = self._d_alpha(i0, i1)
        return self.r ** 2 / 2 * (alpha - np.sin(alpha)) + (self.solve_y(i0) + self.solve_y(i1)) / 2 * (i1 - i0)

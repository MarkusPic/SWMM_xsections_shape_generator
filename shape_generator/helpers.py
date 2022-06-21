import math
import warnings
from abc import ABC, abstractmethod
import numpy as np


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
class _CustomExpr(ABC):
    def __init__(self):
        self.x0 = None
        self.y0 = None

        self.x1 = None
        self.y1 = None

    def __repr__(self):
        return 'Custom Function'

    @abstractmethod
    def expr(self, x):
        pass

    @abstractmethod
    def solve_y(self, i):
        pass

    @abstractmethod
    def solve_x(self, y):
        pass

    @abstractmethod
    def length(self, i0, i1):
        pass

    @abstractmethod
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

    def get_points(self, *args, **kwargs):
        """
        get the point coordinates for the shape-generator

        Args:
            *args:
            **kwargs:

        Returns:
            tuple[list, list]
        """
        pass

####################################################################################################################
class _Linear(_CustomExpr, ABC):
    """abstract linear function"""
    def __init__(self):
        _CustomExpr.__init__(self)

    def set_points(self, start, end):
        """
        set the slope by giving the start and end point

        Args:
            start (tuple[float, float]): start point of the linear function
            end (tuple[float, float]): end point of the linear function
        """
        self.set_start_point(start)
        self.set_end_point(end)

    @classmethod
    def from_points(cls, start, end):
        """
        get a linear function by giving the start and end point

        Args:
            start (tuple[float, float]): start point of the linear function
            end (tuple[float, float]): end point of the linear function

        Returns:
            _Linear: linear function
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

    def get_points(self, *args, **kwargs):
        """get the point coordinates for the shape-generator"""
        x0, y0 = self.get_start_point()
        x1, y1 = self.get_end_point()
        return [x0, x1], [y0, y1]


class Slope(_Linear):
    """
    get function/expression of a straight line with a given point which it intercepts

    Args:
        slope (float): slope
        unit (str): point as a set of a x and a y coordinate

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

        _Linear.__init__(self)

    def __repr__(self):
        return f'Slope Function (k={self.slope:0.2f}, zero=[{self.x0:0.2f}, {self.y0:0.2f}])'

    def expr(self, x):
        """get sympy expression"""
        return self.y0 + (x - self.x0) / self.slope

    def solve_y(self, x):
        """get y value"""
        return self.y0 + (x - self.x0) / self.slope

    def solve_x(self, y):
        """get x value"""
        return self.x0 + (y - self.y0) * self.slope

    def length(self, i0, i1):
        """get shape length between two values"""
        return math.sqrt((self.solve_y(i0) - self.solve_y(i1)) ** 2 + (i0 - i1) ** 2)

    def area(self, x0, x1):
        """get shape area between two values"""
        return (self.solve_y(x0) + self.solve_y(x1)) / 2 * np.abs(x0 - x1)


####################################################################################################################
class Vertical(_Linear):
    """
    function of a vertical line

    for a shape curve this means a constant width
    """
    def __init__(self, y=None):
        """

        Args:
            y (float): y value of the vertical line
        """
        _Linear.__init__(self)
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
class Horizontal(_Linear):
    """
    function of a horizontal line
    """
    def __init__(self, x=None):
        _Linear.__init__(self)
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
class Circle(_CustomExpr):
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

        _CustomExpr.__init__(self)

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
        return np.sqrt(np.maximum(self.r ** 2 - (x - self.x_m) ** 2, 0)) * (-1 if self.clockwise else 1) + self.y_m

    def solve_x(self, y):
        return self.x_m - np.sqrt(self.r ** 2 - (y - self.y_m) ** 2) * (-1 if self.clockwise else 1)

    def length(self, i0, i1):
        return self._d_alpha(i0, i1) * self.r

    def area(self, i0, i1):
        alpha = self._d_alpha(i0, i1)
        return self.r ** 2 / 2 * (alpha - np.sin(alpha)) + (self.solve_y(i0) + self.solve_y(i1)) / 2 * (i1 - i0)

    def get_points(self, step):
        """get the point coordinates for the shape-generator"""
        nx = np.arange(self.x0, self.x1 + step, step).clip(max=self.x1)
        ny = self.solve_y(nx)
        return zip(*ramer_douglas(list(zip(list(nx), list(ny))), dist=step))


####################################################################################################################
class PowerExpr(_CustomExpr):
    def __init__(self, exponent, p=1, x_bottom=0, y_bottom=0):
        """
        power function

        Args:
            exponent (float): exponent
            x_m (float): x axis value of the mid point
            y_m (float): y axis value of the mid point
        """
        self.exponent = float(exponent)
        self.p = p
        # self.x_bottom = float(x_bottom)
        # self.y_bottom = float(y_bottom)

        _CustomExpr.__init__(self)

    def __repr__(self):
        # return f'Power Function (exponent={self.exponent:0.2f}, bottom=[{self.x_bottom:0.2f}, {self.y_bottom:0.2f}])'
        return f'Power Function (exponent={self.exponent:0.2f})'

    def expr(self, x):
        """
        get function/expression of a power function with a given mid point


        Returns:
            sympy.core.expr.Expr: function of the circle
        """
        import sympy as sy
        return sy.Float(self.p) * x ** sy.Float(self.exponent)

    def solve_y(self, x):
        return self.p * x ** (1/self.exponent)

    def solve_x(self, y):
        return (y / self.p) ** self.exponent

    def _first_derivative(self, i):
        return (1/self.exponent) * i ** ((1/self.exponent)-1)

    def length(self, i0, i1):
        from scipy.integrate import quad
        return quad(lambda i: (1 + self._first_derivative(i))**(1/2), i0, i1)[0]

    def area(self, i0, i1):
        from scipy.integrate import quad
        return quad(self.solve_y, i0, i1)[0]

    def get_points(self, step):
        """get the point coordinates for the shape-generator"""
        nx = np.arange(self.x0, self.x1 + step, step).clip(max=self.x1)
        ny = self.solve_y(nx)
        return zip(*ramer_douglas(list(zip(list(nx), list(ny))), dist=step))


####################################################################################################################
# Python package:
# https://github.com/fhirschmann/rdp
#
# ist aber nicht schnell (Faktor 10 langsamer)
#
# https://en.wikipedia.org/wiki/Ramer–Douglas–Peucker_algorithm
#
# Code von hier:
# https://stackoverflow.com/questions/2573997/reduce-number-of-points-in-line#10934629

def _vec2d_dist(p1, p2):
    return (p1[0] - p2[0])**2 + (p1[1] - p2[1])**2


def _vec2d_sub(p1, p2):
    return p1[0]-p2[0], p1[1]-p2[1]


def _vec2d_mult(p1, p2):
    return p1[0]*p2[0] + p1[1]*p2[1]


def ramer_douglas(line, dist):
    """Does Ramer-Douglas-Peucker simplification of a curve with `dist`
    threshold.

    `line` is a list-of-tuples, where each tuple is a 2D coordinate

    Usage is like so:

    >>> myline = [(0.0, 0.0), (1.0, 2.0), (2.0, 1.0)]
    >>> simplified = ramer_douglas(myline, dist = 1.0)
    """

    if len(line) < 3:
        return line

    begin = line[0]
    end = line[-1] if line[0] != line[-1] else line[-2]

    distSq = []
    for curr in line[1:-1]:
        tmp = (
            _vec2d_dist(begin, curr) - _vec2d_mult(_vec2d_sub(end, begin), _vec2d_sub(curr, begin)) ** 2 / _vec2d_dist(begin, end))
        distSq.append(tmp)

    maxdist = max(distSq)
    if maxdist < dist ** 2:
        return [begin, end]
    # print(maxdist)
    pos = distSq.index(maxdist)
    return (ramer_douglas(line[:pos + 2], dist) +
            ramer_douglas(line[pos + 1:], dist)[1:])

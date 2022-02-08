import warnings

import numpy as np

from .helpers import (Circle, _CustomExpr, Slope, Horizontal, get_intersection_point,
                      Vertical, ramer_douglas)

g = 9.81  # m/s² Erdbeschleunigung
ny = 1.31e-6  # m²/s bei 10°C von Wasser

# number of points to describe the shape of the cross-section
# 100 is the limit of points which can be used as a SWMM shape
MAX_NUMBER_POINTS = 100


########################################################################################################################
########################################################################################################################
class CrossSection:
    """main class
    A Class that should help to generate custom cross-section shapes for the SWMM software.

    Attributes:
        label (str): name/label/number of the cross-section
        shape (list): descriptions of the cross-section as commands in a list
        accuracy (int): number of decimal points to use for the export
        unit (str): unit of entered values
        double (bool): if the cross-section two separate cross-sections
        points (list): points of the final cross-section
        simplify (bool): if the ramer-douglas algorithm should be used to simplify curve
    """

    def __init__(self, label, height=None, width=None, unit=None, double=False, accuracy=3, simplify=True):
        """Initialise the cross-section class

        Args:
            label (str): name/label/number of the cross-section
            height (float): absolute height of the CS
            width (Optional[float]): absolute width of the CS (optional) can be estimated
            unit (Optional[str]): enter unit to add the unit in the plots
            simplify (bool): if the ramer-douglas algorithm should be used to simplify curve
        """
        self.label = label

        self._height = height
        self._width = width
        self.unit = unit
        # _______________________________
        self.shape = list()
        self._shape_description = None  # functions to describe the cross-section
        # _______________________________
        self.accuracy = accuracy
        self.double = double
        self.simplify = simplify

        # _______________________________
        # Profile data
        self.points = list()

        # _______________________________
        # calculate stationary flow
        self._area_v = None
        self._r_hyd_v = None
        self._l_u_v = None
        self._r_hyd_v = None
        self._v_v = None

    def __repr__(self):
        return f'CrossSection({self})'

    def __str__(self):
        return f'{self.title_string()}'

    @property
    def identifier(self):
        return f'{self.label}'

    def name_string(self):
        return f'{self.label}'

    def title_string(self):
        return self.label

    @property
    def height(self):
        """
        absolute height of the CS

        Returns:
            float: absolute height of the CS
        """
        return self._height

    @property
    def width(self):
        """
        absolute width of the CS

        Returns:
            float: absolute width of the CS
        """
        return self._width

    def get_width(self):
        """
        get absolute width of cross-section

        Returns:
            float: width of cross-section
        """
        if not self.get_points():
            return None
        else:
            return max(self.points[1]) * 2

    def get_height(self):
        """
        get absolute height of cross-section

        Returns:
            float: height of cross-section
        """
        if not self.get_points():
            return None
        else:
            return max(self.points[0])

    def set_double_cross_section(self):
        """
        make the cross-section as a double section (=Doppelprofil)
        """
        self.double = True

    def _reset_shape(self):
        self.points = list()
        self._shape_description = None

    def add(self, x_or_expr, y=None):
        """
        add part of cross-section

        can be a:

        - function/expression
        - point (x,y) coordinates
        - boundary condition (x or y) of a surrounding function = only x or y is given and the other is :obj:`None`
        - slope (x=slope, y=unit of slope)


        Args:
            x_or_expr (Optional[float , None, CustomExpr]):

                - :obj:`float` : x coordinate or x-axis boundary or slope if any str keyword is used in argument ``y``
                - :obj:`CustomExpr` : Expression/function for the cross-section part,
                  i.e.:  :obj:`shape_generator.Circle`, :obj:`Slope`, :obj:`Vertical`, :obj:`Horizontal`
                - :obj:`None` : if a y-axis boundary is given

            y (Optional[float,str]): y coordinate of unit of slope

                - :obj:`float` : x coordinate or x-axis boundary
                - :obj:`None` : if a x-axis boundary is given or an expression in ``x_or_expr``
                - :obj:`str` : argument x is a slope

                    - ``slope`` : ready to use slope 1 / :math:`\\Delta` y
                    - ``°slope`` : slope as an angle in degree (°)
                    - ``%slope`` : slope in percentage (%)

        """
        self._reset_shape()

        if isinstance(x_or_expr, _CustomExpr):
            self.shape.append(x_or_expr)

        else:
            if x_or_expr is not None:
                x = float(x_or_expr)
            else:
                x = x_or_expr

            if isinstance(y, str) and 'slope' in y:
                if x == 0:
                    self.shape.append(Horizontal())
                else:
                    unit = y.replace('slope', '')
                    self.shape.append(Slope(x, unit=unit))

            else:
                if y is not None:
                    y = float(y)

                self.shape.append((x, y))

    @property
    def shape_description(self):
        """
        list of functions to describe the cross-section shape

        Returns:
            list: description of the cross-section shape
        """
        if self._shape_description is None:
            # result list
            expr_list = []

            # boundary condition
            point_final = (self.height, 0)

            is_next_final_point = False
            for i, shape_i in enumerate(self.shape):
                # _________________________________
                if i == 0:
                    point_prev = (0, 0)
                else:
                    point_prev = expr_list[-1].get_end_point()

                # _________________________________
                # boundary condition
                if (i + 1) == len(self.shape):
                    is_next_final_point = True
                    shape_next = point_final
                else:
                    shape_next = self.shape[i + 1]

                # ____________________________________________________________
                if isinstance(shape_i, tuple):
                    # shape_i ist ein tuple aus x und y Koordinate (x,y)
                    # nächster Punkt in der Shape
                    point_i = shape_i
                    expr_i = None

                    if (point_i[0] is None) or (point_i[1] is None):
                        # this part is only used as boundary condition
                        if is_next_final_point:
                            expr_i = Slope.from_points(point_prev, point_final)
                            expr_list.append(expr_i)
                        continue

                    if (point_prev[1] is not None) and (point_prev != point_i):
                        expr_i = Slope.from_points(point_prev, point_i)
                        expr_list.append(expr_i)

                    if is_next_final_point and (point_i != point_final):
                        expr_i = Slope.from_points(point_i, point_final)
                        expr_list.append(expr_i)

                    del point_i, expr_i

                # ____________________________________________________________
                elif isinstance(shape_i, _CustomExpr):
                    # shape_i ist eine Funktion zur Beschreibung der Shape
                    expr_i = shape_i  # CustomExpr = Horizontal, Vertical, Circle, Slope

                    expr_i.set_start_point(point_prev)

                    if isinstance(shape_next, tuple):
                        expr_i.set_end_point(shape_next)

                    elif isinstance(shape_next, _CustomExpr):
                        expr_i.set_end_point(get_intersection_point(expr_i, expr2=shape_next,
                                                                    x_from=point_prev[0], x_to=self.height))

                    else:
                        raise NotImplementedError('Unknown Input in shape')

                    expr_list.append(expr_i)

                    del expr_i

            # ____________________________________________________________
            self._shape_description = expr_list

        return self._shape_description

    def iter_shape_description(self):
        for shape in self.shape_description:  # type: _CustomExpr
            yield shape.x0, shape.x1, shape

    def get_points(self):
        """create absolute point coordinates and write it into :py:attr:`~points`

        To create a :obj:`list[tuple]` of all the points to describe the cross-section.
        This function replaces the Expressions given in :py:attr:`~add` to points with x and y coordinates
        and writes them into the :py:attr:`~points` attribute.

        Returns:
            list[list[float,float]]: absolute point coordinates
        """
        if not self.points:
            step = 10 ** (-self.accuracy) * self.height
            x = list()
            y = list()

            for shape in self.shape_description:
                x_i, y_i = shape.get_points(step)
                # drop duplicate point in intersections
                if x and (x_i[0] == x[-1]) and (y_i[0] == y[-1]):
                    x += x_i[1:]
                    y += y_i[1:]
                else:
                    x += x_i
                    y += y_i

            # ----------
            if self.simplify:
                x, y = zip(*ramer_douglas(list(zip(x, y)), dist=step))
            else:
                x_, y_ = zip(*ramer_douglas(list(zip(x, y)), dist=step))
                if len(x) != len(x_):
                    print(f'Possible reduction with keyword `simplify=True´: from {len(x)} to {len(x_)} points '
                          f'(with Ramer-Douglas, distance={step})')

            # ----------
            if len(x) > MAX_NUMBER_POINTS:

                if not self.simplify and (len(x_) <= MAX_NUMBER_POINTS):
                    warnings.warn(f'to many points (n={len(x)}) -> set simplify=True (n={len(x_)})')
                    self.simplify = True
                    return self.get_points()

                # ----------
                # number of expressions used in shape
                warnings.warn(f'to many points (n={len(x)}) -> reduce accuracy')
                num_arcs = sum([isinstance(i, Circle) for i in self.shape_description])

                if not num_arcs:
                    raise UserWarning('No arcs but to many points -> UNKNOWN')

                # number of fixed points in shape
                num_points = (len(self.shape_description) - num_arcs) * 2

                # calculate the net height of the circle functions.
                function_steps = {i: s[1] - s[0] for i, s in enumerate(self.iter_shape_description()) if
                                  isinstance(self.shape_description[i], Circle)}
                # step size used to discretise the expressions
                step2 = sum(function_steps.values()) / (MAX_NUMBER_POINTS - num_points)

                self.accuracy = np.log10(step2 / self.height) * -1

                return self.get_points()

            # -------------------------
            # prevent duplicate x values (raises SWMM error)
            if len(x[1:-1]) != len(set(x[1:-1])):
                x = list(x)
                for i in range(1, len(x)-1):
                    if (x[i] != 0) and (x[i] == x[i-1]):
                        x[i] += step

            # -------------------------
            self.points = x, y

        return self.points

    def profile_axis(self, ax, relative=False, half=False, fill=False, marker='.', ls='-', **kwargs):
        """
        plot the shape curve on an axes

        Args:
            ax (matplotlib.pyplot.Axes): plot axes
            relative (bool): if the plot should be in relative size
            half (bool): if ony the half curve should be plotted
            fill (bool): if the curve should be fill on the inside
            marker (str): marker of the curve
            ls (str): line-style of the curve
            **kwargs: plot keyword arguments

        Returns:
            matplotlib.pyplot.Axes: plot axes
        """
        x, y = self.get_points()
        hi = np.array(x)
        wi = np.array(y)

        height = self.get_height()

        if relative:
            hi /= height
            wi /= height

        if not half:
            hi = np.append(hi, hi[::-1])
            wi = np.append(wi, wi[::-1]*-1)

        # -------------------------
        ax.plot(wi, hi, marker=marker, ls=ls, zorder=1000000, clip_on=False, **kwargs)
        if fill:
            ax.fill(wi, hi)

        return ax

    def profile_figure(self, relative=False, half=False, fill=False, **kwargs):
        """
        create a plot of the cross-section

        Args:
            relative (bool): if the plot should be in relative size
            half (bool): if ony the half curve should be plotted
            fill (bool): if the curve should be fill on the inside
            **kwargs: plot keyword arguments

        Returns:
            matplotlib.pyplot.Figure: figure of the plot
        """
        import matplotlib.pyplot as plt

        def ceil_base(i, base):
            return base * np.ceil(float(i) / base)

        # -------------------------
        fig, ax = plt.subplots()

        ax = self.profile_axis(ax, relative=relative, half=half, fill=fill, **kwargs)
        h = self.get_height()
        w = self.get_width() / 2
        # -------------------------
        title = self.title_string()

        if relative:
            base = 0.1
            w /= h
            h = 1
            # -------------------------
            ax.set_ylabel('rel H')
            ax.set_xlabel('$^B/_H$')

        else:
            base = 10 ** np.floor(np.log10(w))

            if base % 1 == 0:
                base = int(base)

            # -------------------------
            title = f'{title}\n{h}x{ceil_base(w * 2, base/2)}'
            if self.unit is not None:
                title += self.unit

        ax.set_title(title)
        xlim = ceil_base(w, base)
        ylim = ceil_base(h, base)

        # title += f'\nRaster: {base}'
        ax.text(xlim - base * .5, 0, str(base), size='xx-small', ha='center', va='top')
        ax.text(xlim, base * .5, str(base), size='xx-small', va='center', ha='left', rotation=90)

        # -------------------------
        if half:
            xlim_left = 0
        else:

            xlim_left = -xlim

        # -------------------------
        ax.set_aspect('equal', 'box')
        ax.set_xticks(np.arange(xlim_left, xlim, base), minor=False)
        if base / 2 != 0:
            ax.set_xticks(np.arange(xlim_left, xlim, base / 2), minor=True)

        ax.set_yticks(np.arange(0, ylim, base), minor=False)
        if base / 2 != 0:
            ax.set_yticks(np.arange(0, ylim, base / 2), minor=True)

        ax.tick_params(which='both', length=0, width=0, labelbottom=False, labeltop=False, labelleft=False,
                       labelright=False, bottom=False, top=False, left=False, right=False)

        ax.set_xlim(xlim_left, xlim)
        ax.set_ylim(0, ylim)
        ax.grid(True)
        ax.set_axisbelow(True)

        # ------------------
        fig.tight_layout()
        return fig

    ####################################################################################################################
    # testing new functions
    ####################################################################################################################
    def b_w_t(self, hi):
        """
        width of the cross-section at a certain height

        (Wasseroberflächenbreite im teilgefüllten Querschnitt)

        Args:
            hi (float | numpy.ndarray): a certain height

        Returns:
            float | numpy.ndarray: width at the certain height
        """
        if isinstance(hi, np.ndarray):
            w = np.array([np.NaN] * hi.size)
            # w = hi.copy()

            for i, (lower, upper, f) in enumerate(self.iter_shape_description()):
                b = (hi >= lower) & (hi <= upper)
                w[b] = f.solve_y(hi[b])

            return w * 2
        else:
            for lower, upper, f in self.iter_shape_description():
                if lower <= hi <= upper:
                    return f.solve_y(hi) * 2

    def l_u_t(self, hi):
        """
        wetted perimeter in the partial-filled cross-section at a certain water level height

        (benetzter Umfang im teilgefüllten Querschnitt)

        Args:
            hi (float | numpy.ndarray): a certain height

        Returns:
            float | numpy.ndarray: wetted perimeter at the certain height
        """
        if isinstance(hi, np.ndarray):
            l = np.array([0.] * hi.size)

            for i, (lower, upper, f) in enumerate(self.iter_shape_description()):
                b = hi > upper
                l[b] += f.length(lower, upper)

                b = (hi >= lower) & (hi <= upper)
                l[b] += f.length(lower, hi[b])

        else:
            l = 0
            for lower, upper, f in self.iter_shape_description():
                if hi > upper:
                    l += f.length(lower, upper)
                elif lower <= hi <= upper:
                    l += f.length(lower, hi)
                    break
                else:
                    break

        return l * 2

    @property
    def l_u_v(self):
        """
        wetted perimeter of the full-filled cross-section

        (benetzter Umfang im vollgefüllten Querschnitt)

        Returns:
            float | numpy.ndarray: wetted perimeter
        """
        if self._l_u_v is None:
            self._l_u_v = self.l_u_t(self.height)
        return self._l_u_v

    def area_t(self, hi):
        """
        flow area in the partial-filled cross-section at a certain water level height

        (Fließquerschnitt im teilgefüllten Querschnitt)

        Args:
            hi (float | numpy.ndarray): a certain height

        Returns:
            float | numpy.ndarray: flow area at the certain height
        """
        if isinstance(hi, np.ndarray):
            a = np.array([0.] * hi.size)

            for i, (lower, upper, f) in enumerate(self.iter_shape_description()):
                b = hi > upper
                a[b] += f.area(lower, upper)

                b = (hi >= lower) & (hi <= upper)
                a[b] += f.area(lower, hi[b])

        else:
            a = 0
            for lower, upper, f in self.iter_shape_description():
                if hi > upper:
                    a += f.area(lower, upper)
                elif lower <= hi <= upper:
                    a += f.area(lower, hi)
                    break
                else:
                    break

        return a * 2

    @property
    def area_v(self):
        """
        flow area of the full-filled cross-section

        (Fließquerschnitt im vollgefüllten Querschnitt)

        Returns:
            float | numpy.ndarray: flow area
        """
        if self._area_v is None:
            self._area_v = self.area_t(self.height)
        return self._area_v

    def r_hyd_t(self, hi):
        """
        hydraulic radius in the partial-filled cross-section at a certain water level height

        (hydraulischer Radius im teilgefüllten Querschnitt)

        Args:
            hi (float | numpy.ndarray): a certain height

        Returns:
            float | numpy.ndarray: hydraulic radius at the certain height
        """
        return self.area_t(hi) / self.l_u_t(hi)

    @property
    def r_hyd_v(self):
        """
        hydraulic radius of the full-filled cross-section

        (hydraulischer Radius im vollgefüllten Querschnitt)

        Returns:
            float | numpy.ndarray: hydraulic radius
        """
        if self._r_hyd_v is None:
            self._r_hyd_v = self.area_v / self.l_u_v
        return self._r_hyd_v

    ####################################################################################################################
    def velocity_v(self, slope, k):
        """
        velocity in a full-filled channel

        Args:
            slope (float): absolute slope in m/m
            k (float): Betriebliche Rauigkeit in mm

        Returns:
            float: full-filling velocity in m/s

        References:
            DWA-A 110 Section 4.1.1 Vollfüllung
        """
        if self._v_v is None:
            self._v_v = dict()

        if k not in self._v_v:
            self._v_v[k] = dict()

        if slope not in self._v_v[k]:
            self._v_v[k][slope] = None

        if self._v_v[k][slope] is None:
            r_hyd = self.r_hyd_v / 1000  # from mm to m
            J = slope  # / 1000
            k = k / 1000  # from mm to m
            self._v_v[k][slope] = (
                    -2 * np.log10(2.51 * ny / (4 * r_hyd * np.sqrt(2 * g * J)) + k / (14.84 * r_hyd)) * np.sqrt(
                2 * g * 4 * r_hyd * J))

        return self._v_v[k][slope]

    def velocity_t(self, hi, slope, k):
        """
        velocity in a partial-filled channel at a certain water level height

        Args:
            hi (float): water level = height in mm
            slope (float): absolute slope in m/m
            k (float): Betriebliche Rauigkeit in mm

        Returns:
            float: velocity in m/s

        References:
            DWA-A 110 Section 4.1.2 Teilfüllung
        """
        return (self.r_hyd_t(hi) / self.r_hyd_v) ** 0.625 * self.velocity_v(slope, k)

    def flow_t(self, hi, slope, k):
        """
        flow in a partial-filled channel at a certain water level height

        Args:
            hi (float): water level = height in mm
            slope (float): absolute slope in m/m
            k (float): Betriebliche Rauigkeit in mm

        Returns:
            float: flow rate in L/s
        """
        return self.velocity_t(hi, slope, k) * self.area_t(hi) * 1.0e-6 * 1000

    def flow_v(self, slope, k):
        """
        flow in a full-filled channel

        Args:
            slope (float): absolute slope in m/m
            k (float): Betriebliche Rauigkeit in mm

        Returns:
            float: full-filled flow rate in L/s
        """
        return self.velocity_v(slope, k) * self.area_v * 1.0e-6 * 1000

    ####################################################################################################################
    def h_t(self, Q_t, slope, k):
        """
        get the height of the water level based on the known flow

        Args:
            Q_t (float): flow in L/s
            slope (float): absolute slope in m/m
            k (float): Betriebliche Rauhigkeit in mm

        Returns:
            float: height of the water level
        """
        # hi = '?'
        # self.flow_t(hi, slope, k)
        # self.flow_v(slope, k)
        from scipy.optimize import minimize_scalar
        res = minimize_scalar(lambda hi: abs(Q_t - self.flow_t(hi, slope, k)), bounds=(0, self.height),
                              method='bounded')
        return res.x

    ####################################################################################################################
    # swmm_api functions
    ####################################################################################################################
    def to_curve(self):
        """create a SWMM curve object with the data of the swmm-shape_generator-CrossSection"""
        from swmm_api.input_file.sections import Curve
        x, y = self.get_points()
        # [1:-1] without first and last point / not needed in swmm
        height = np.array(x[1:-1]) / self.height
        area = np.array(y[1:-1]) / self.height * 2
        return Curve(Name=self.identifier, Type=Curve.TYPES.SHAPE,
                     points=[[float(h), float(a)] for h, a in zip(height, area)])

    @classmethod
    def from_curve(cls, curve, height=100, *args, **kwargs):
        """
        create an object with the data of the swmm curve data as relative coordinates

        Args:
            curve (swmm_api.input_file.sections.Curve): Curve object of the CURVES section in the inp-data file
            height (float): absolute height of the CS
            *args: arguments, see :attr:`CrossSection.__init__`
            **kwargs: keyword arguments, see :attr:`CrossSection.__init__`

        Returns:
            CrossSection: of the shape coordinates
        """
        cross_section = CrossSection(curve.Name, height=height, *args, **kwargs)
        for x, y in curve.points:
            cross_section.add(x * height, y * height / 2)
        return cross_section

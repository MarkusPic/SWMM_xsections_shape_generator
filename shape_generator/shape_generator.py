import warnings
from os import path

import matplotlib.pyplot as plt
import sympy as sy
from numpy import array, arange, ndarray, ceil, log10, floor, NaN, append

from .curve_simplification import ramer_douglas
from .helpers import Circle, x, CustomExpr, Slope, Horizontal, sqrt

g = 9.81  # m/s^2 Erdbeschleunigung
ny = 1.31e-6  # m^2/s bei 10°C von Wasser


########################################################################################################################
########################################################################################################################
class CrossSection:
    """main class
    A Class that should help to generate custom cross section shapes for the SWMM software.

    Attributes:
        label (str): name/label/number of the cross section
        description (Optional[str]): optional description of the cross section
        shape (list): descriptions of the cross section as commands in a list
        accuracy (int): number of decimal points to use for the export
        working_directory (str): directory where the files get saved
        unit (str): unit of entered values
        double (bool): if the cross section two separate cross sections
    """

    def __init__(self, label, description=None, height=None, width=None, working_directory='', unit=None):
        """Initialise the cross section class

        Args:
            label (str): main name/label/number of the cross section
            description (Optional[str]): optional description of the cross section
            height (float): absolute height of the CS
            width (Optional[float]): absolute width of the CS (optional) can be estimated
            working_directory (str): directory where the files are saved
            unit (Optional[str]): enter unit to add the unit in the plots
        """
        self.label = label
        self.description = ''
        if description is not None:
            self.description = str(description).strip()

        self._height = height
        self._width = width
        self.shape = list()
        self._shape_description = None  # functions to describe the cross section
        self.accuracy = 3
        self.working_directory = working_directory
        self.unit = unit
        self.double = False

        # _______________________________
        # Profile data
        self.points = list()

        # _______________________________
        # print('_' * 30)
        # print(self)

        # _______________________________
        # calculate stationary flow
        self._area_v = None
        self._r_hyd_v = None
        self._l_u_v = None
        self._r_hyd_v = None
        self._v_v = None

        # _______________________________
        # number of points to describe the shape of the cross section
        # 100 is the limit of points which can be used as a SWMM shape
        self.max_number_points = 100

    def __repr__(self):
        return str(self)

    def __str__(self):
        return '{}:  {}'.format(self.label, self.description)

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
        get absolute width of cross section

        Returns:
            float: width of cross section
        """
        if not self.get_points():
            return None
        else:
            return max(self.points[1]) * 2

    @property
    def out_filename(self):
        """
        filename of the figure/text-file to be created

        Returns:
            str: filename

        """
        return path.join(self.working_directory, str(self.label))

    def _reset_shape(self):
        self.points = list()
        self._shape_description = None

    def add(self, x_or_expr, y=None):
        """
        add part of cross section

        can be a:

        - function/expression
        - point (x,y) coordinates
        - boundary condition (x or y) of a surrounding function = only x or y is given and the other is :obj:`None`
        - slope (x=slope, y=unit of slope)


        Args:
            x_or_expr (Optional[float , None, CustomExpr]):

                - :obj:`float` : x coordinate or x-axis boundary or slope if any str keyword is used in argument ``y``
                - :obj:`CustomExpr` : Expression/function for the cross section part,
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

        if isinstance(x_or_expr, CustomExpr):
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

    def set_double_cross_section(self):
        """
        make the cross section as a double section (=Doppelprofil)
        """
        self.double = True

    @property
    def shape_description(self):
        """
        list of functions to describe the cross section shape

        Returns:
            list: description of the cross section shape
        """
        if self._shape_description is None:
            # result list
            function = list()

            def add_slope_to_function(point0, point1):
                start = point0[0]
                end = point1[0]
                yi = Slope.from_points(point0, point1)
                function.append((start, end, yi))

            # boundary condition
            last_point = (0, 0)
            final_point = (self.height, 0)

            for i, shape_i in enumerate(self.shape):

                # _________________________________
                # boundary condition
                if (i + 1) == len(self.shape):
                    shape_next = final_point
                else:
                    shape_next = self.shape[i + 1]

                # _________________________________
                # if isinstance(shape_i, tuple) and shape_i[1] == 'slope':
                #     shape_i = Slope(shape_i[0])
                #     shape_i.set_start_point(last_point)

                # ____________________________________________________________
                if isinstance(shape_i, tuple):

                    if (shape_i[0] is None) or (shape_i[1] is None):
                        # this part is only used as boundary condition
                        if shape_next == final_point:
                            start = last_point[0]
                            end = shape_next[0]
                            yi = Slope.from_points(last_point, shape_next)
                            function.append((start, end, yi))

                        continue

                    if last_point[1] is not None:
                        start = last_point[0]
                        end = shape_i[0]
                        yi = Slope.from_points(last_point, shape_i)
                        function.append((start, end, yi))

                    if shape_next == final_point:
                        start = shape_i[0]
                        end = shape_next[0]
                        yi = Slope.from_points(shape_i, shape_next)
                        function.append((start, end, yi))

                    # ________________________________

                    last_point = (end, shape_i[1])

                # ____________________________________________________________
                elif isinstance(shape_i, CustomExpr):
                    yi = shape_i

                    if isinstance(yi, Slope) and yi.x0 is None:
                        yi.set_start_point(last_point)

                    start = last_point[0]

                    if isinstance(yi, Horizontal):
                        if isinstance(shape_next, tuple):
                            yi.set_points(last_point, shape_next)

                        elif isinstance(shape_next, CustomExpr):
                            warnings.warn('must be implemented', FutureWarning)

                    else:
                        if isinstance(shape_next, tuple):
                            end = shape_next[0]

                            if end is None and shape_next[1] is not None:
                                end = sy.solve(yi.expr() - shape_next[1], x)[0]

                        elif isinstance(shape_next, CustomExpr):
                            res = sy.solve(yi.expr() - shape_next.expr(), x)
                            if len(res) == 0:
                                from scipy.optimize import minimize_scalar
                                end = minimize_scalar(lambda j: float((yi.expr() - shape_next.expr()).subs(x, j)),
                                                      bounds=(start, self.height), method='bounded').x

                            elif len(res) == 1:
                                end = float(res[0])
                            else:
                                # multiple results
                                # TODO: how to handle it
                                end = float(res[0])

                        else:
                            raise NotImplementedError('Unknown Input in shape')

                        end = float(end)

                        if start == end:
                            warnings.warn('unused part of the shape detected. Ignoring this part.')
                            continue

                    function.append((start, end, yi))

                    # ____________________________
                    if isinstance(shape_next, tuple) and shape_next[1] is not None:
                        last_point = (end, shape_next[1])
                    else:
                        last_point = (end, float(yi.solve(end)))

            # ____________________________________________________________
            self._shape_description = function

        return self._shape_description

    def _get_points_legacy(self):
        """create absolute point coordinates and write it into :py:attr:`~df_abs`

        To create a :obj:`pandas.DataFrame` of all the points to describe the cross section.
        This function replaces the Expressions given in :py:attr:`~add` to points with x and y coordinates
        and writes them into the :py:attr:`~df_abs` attribute.

        Returns:
            pandas.DataFrame: absolute point coordinates
        """
        # number of expressions used in shape
        warnings.warn('get_points | legacy mode')
        num_functions = sum([isinstance(i[2], Circle) for i in self.shape_description])

        step = None
        # if functions are used in shape

        if num_functions:
            # number of fixed points in shape
            num_points = (len(self.shape_description) - num_functions) * 2

            # calculate the net height of the circle functions.
            function_steps = {i: s[1] - s[0] for i, s in enumerate(self.shape_description) if
                              isinstance(self.shape_description[i][2], Circle)}
            # step size used to discretise the expressions
            step = sum(function_steps.values()) / (self.max_number_points - num_points)

            min_step = 1 * 10 ** (-self.accuracy) * self.height
            if step < min_step:
                step = min_step

        x = list()
        y = list()

        for start, end, f in self.shape_description:
            if isinstance(f, Circle):
                # this_step = (end - start) / floor((end - start) / step)
                # print(step, ' vs ', this_step)
                nx = arange(start, end + step, step).clip(max=end)
                ny = f.solve(nx)
                x += list(nx)
                y += list(ny)
            elif isinstance(f, Horizontal):
                x0, y0 = f.start_point()
                x1, y1 = f.end_point()
                x += [x0, x1]
                y += [y0, y1]
            else:
                nx = array([start, end])
                x += list(nx)
                y += list(f.solve(nx))

        return x, y

    def get_points(self):
        """create absolute point coordinates and write it into :py:attr:`~points`

        To create a :obj:`list[tuple]` of all the points to describe the cross section.
        This function replaces the Expressions given in :py:attr:`~add` to points with x and y coordinates
        and writes them into the :py:attr:`~points` attribute.

        Returns:
            list[list[float,float]]: absolute point coordinates
        """
        if not self.points:
            step = 10 ** (-self.accuracy) * self.height
            # if functions are used in shape

            x = list()
            y = list()

            for start, end, f in self.shape_description:
                if isinstance(f, Circle):
                    nx = arange(start, end + step, step).clip(max=end)
                    ny = f.solve(nx)
                    x += list(nx)
                    y += list(ny)
                elif isinstance(f, Horizontal):
                    x0, y0 = f.start_point()
                    x1, y1 = f.end_point()
                    x += [x0, x1]
                    y += [y0, y1]
                else:
                    nx = array([start, end])
                    x += list(nx)
                    y += list(f.solve(nx))
            x, y = zip(*ramer_douglas(list(zip(x, y)), dist=step))

            if len(x) > self.max_number_points:
                x, y = self._get_points_legacy()

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

    def _check_points(self):
        """
        remove errors from the point cloud, ie.:
            - remove duplicates,
            - (if specified) remove points which overlap the overall cross section width and
            - other errors...
        """
        df = self.df_rel  # (self.df_abs / self.height).copy()
        df = df.round(self.accuracy)
        df = df.drop_duplicates()

        if self.width is not None and any(df['y'] > self.width / 2):
            df['y'] = df['y'].clip_upper(self.width / 2)
            warnings.warn('had to clip the width')

        # print((arctan(df['x'].diff() / df['y'].diff())/ pi * 180).head(10))
        df = df.dropna()

        if self.double:
            df['y'] *= 2

        # delete errors
        df = df[df['x'].expanding().max() == df['x']].copy()

        # # change x duplicates
        # dupls = df['x'].duplicated(keep=False)
        # if dupls.any():
        #     nx = df['x'][dupls]
        #
        #     def raise_values(s):
        #         return s + Series(index=s.index, data=range(len(s.index))) * 10 ** (-self.accuracy)
        #
        #     nx = nx.groupby(nx).apply(raise_values)
        #     df.loc[nx.index, 'x'] = nx

        self._df_abs = (df * self.height).copy()

    def profile_axis(self, ax, relative=False, half=False, fill=False, marker='.', ls='-', **kwargs):
        x, y = self.get_points()
        hi = array(x)
        wi = array(y)

        w = wi.max()
        h = hi.max()

        if relative:
            hi /= h
            wi /= h

        if not half:
            hi = append(hi, hi[::-1])
            wi = append(wi, wi[::-1]*-1)

        # -------------------------
        ax.plot(wi, hi, marker=marker, ls=ls, zorder=1000000, clip_on=False, **kwargs)
        if fill:
            ax.fill(wi, hi)

        return ax, (h, w)

    def profile_figure(self, relative=False, half=False, fill=False, **kwargs) -> plt.Figure:
        """create a plot of the cross section"""
        def custom_round(x_, base):
            return base * ceil(float(x_) / base)

        # -------------------------
        fig, ax = plt.subplots()

        ax, (h, w) = self.profile_axis(ax, relative=relative, half=half, fill=fill, **kwargs)
        # -------------------------
        if relative:
            xlim = 1
            ylim = 1
            base = 0.2

            # -------------------------
            ax.set_ylabel('rel H')
            ax.set_xlabel('B/H')
            ax.set_title('{}: {}'.format(self.label, self.description))

        else:
            base = 10 ** floor(log10(w))
            xlim = custom_round(w, base)
            ylim = custom_round(h, base)

            # -------------------------
            n = self.label
            if self.label != self.description:
                n += ': {}'.format(self.description)

            ax.set_title('{}\n{:0.0f}x{:0.0f}'.format(n, h, custom_round(w * 2, base/2)) +
                         (self.unit if self.unit is not None else ''))

        # -------------------------
        if half:
            xlim_left = 0
        else:

            xlim_left = -xlim

        # -------------------------
        # ax.legend().remove()
        ax.set_aspect('equal', 'box')
        ax.set_xticks(arange(xlim_left, xlim, base), minor=False)
        if base / 2 != 0:
            ax.set_xticks(arange(xlim_left, xlim, base / 2), minor=True)

        ax.set_yticks(arange(0, ylim, base), minor=False)
        if base / 2 != 0:
            ax.set_yticks(arange(0, ylim, base / 2), minor=True)

        # ax.set_axis_off()
        # ax.set_frame_on(False)
        # ax.axis()
        ax.tick_params(which='both', length=0, width=0, labelbottom=False, labeltop=False, labelleft=False,
                       labelright=False, bottom=False, top=False, left=False, right=False)

        ax.set_xlim(xlim_left, xlim)
        ax.set_ylim(0, ylim)
        ax.grid(True)
        # ax.grid(True, which='minor', linestyle=':', linewidth=0.5)
        ax.set_xlabel(None)
        ax.set_axisbelow(True)

        # ------------------
        fig.tight_layout()
        return fig

    ####################################################################################################################
    # testing new functions
    ####################################################################################################################
    def b_w_t(self, hi):
        """
        width of the cross section at a certain height
        (Wasseroberflächenbreite im teilgefüllten Querschnitt)

        Args:
            hi (float | numpy.ndarray): a certain height

        Returns:
            float | numpy.ndarray: width at the certain height
        """
        if isinstance(hi, ndarray):
            w = array([NaN] * hi.size)
            # w = hi.copy()

            for i, (lower, upper, f) in enumerate(self.shape_description):
                b = (hi >= lower) & (hi <= upper)
                w[b] = f.solve(hi[b])

            return w * 2
        else:
            for lower, upper, f in self.shape_description:
                if lower <= hi <= upper:
                    return f.solve(hi) * 2

    def l_u_t(self, hi):
        """
        wetted perimeter in the partially filled cross section at a certain water level height
        (benetzter Umfang im teilgefüllten Querschnitt)

        Args:
            hi (float | numpy.ndarray): a certain height

        Returns:
            float | numpy.ndarray: wetted perimeter at the certain height
        """
        if isinstance(hi, ndarray):
            l = array([0.] * hi.size)

            for i, (lower, upper, f) in enumerate(self.shape_description):
                b = hi > upper
                l[b] += f.length(lower, upper)

                b = (hi >= lower) & (hi <= upper)
                l[b] += f.length(lower, hi[b])

        else:
            l = 0
            for lower, upper, f in self.shape_description:
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
        wetted perimeter of the full filled cross section
        (benetzter Umfang im vollgefüllten Querschnitt)

        Returns:
            float | numpy.ndarray: wetted perimeter
        """
        if self._l_u_v is None:
            self._l_u_v = self.l_u_t(self.height)
        return self._l_u_v

    def area_t(self, hi):
        """
        flow area in the partially filled cross section at a certain water level height
        (Fließquerschnitt im teilgefüllten Querschnitt)

        Args:
            hi (float | numpy.ndarray): a certain height

        Returns:
            float | numpy.ndarray: flow area at the certain height
        """
        if isinstance(hi, ndarray):
            a = array([0.] * hi.size)

            for i, (lower, upper, f) in enumerate(self.shape_description):
                b = hi > upper
                a[b] += f.area(lower, upper)

                b = (hi >= lower) & (hi <= upper)
                a[b] += f.area(lower, hi[b])

        else:
            a = 0
            for lower, upper, f in self.shape_description:
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
        flow area of the full filled cross section
        (Fließquerschnitt im vollgefüllten Querschnitt)

        Returns:
            float | numpy.ndarray: flow area
        """
        if self._area_v is None:
            self._area_v = self.area_t(self.height)
        return self._area_v

    def r_hyd_t(self, hi):
        """
        hydraulic radius in the partially filled cross section at a certain water level height
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
        hydraulic radius of the full filled cross section
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
        calculate velocity in partially filled sewer channel

        Args:
            slope (float): ablosute slope in m/m
            k (float): Betriebliche Rauhigkeit in mm

        Returns:
            float: full filling velocity in m/s

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
                    -2 * log10(2.51 * ny / (4 * r_hyd * sqrt(2 * g * J)) + k / (14.84 * r_hyd)) * sqrt(
                2 * g * 4 * r_hyd * J))

        return self._v_v[k][slope]

    def velocity_t(self, hi, slope, k):
        """
        calculate velocity in partially filled sewer channel

        Args:
            hi (float): water level = height in mm
            slope (float): ablosute slope in m/m
            k (float): Betriebliche Rauhigkeit in mm

        Returns:
            float: velocity in m/s

        References:
            DWA-A 110 Section 4.1.2 Teilfüllung
        """
        return (self.r_hyd_t(hi) / self.r_hyd_v) ** 0.625 * self.velocity_v(slope, k)

    def flow_t(self, hi, slope, k):
        """

        Args:
            hi (float): water level = height in mm
            slope (float): ablosute slope in m/m
            k (float): Betriebliche Rauhigkeit in mm

        Returns:
            float: flow rate in L/s
        """
        return self.velocity_t(hi, slope, k) * self.area_t(hi) * 1.0e-6 * 1000

    def flow_v(self, slope, k):
        """

        Args:
            slope (float): absolute slope in m/m
            k (float): Betriebliche Rauhigkeit in mm

        Returns:
            float: full filling flow rate in L/s
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

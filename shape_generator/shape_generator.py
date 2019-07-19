import warnings
import matplotlib.pyplot as plt
import numpy as np
import sympy as sy
import pandas as pd

from math import radians, cos, ceil, log10, floor
from os import path
from webbrowser import open as open_file
from numpy import NaN
from pandas import isna, notna

from .helpers import deg2slope, channel_end, Circle, x, CustomExpr, Slope, Vertical, Horizontal, sqrt

g = 9.81  # m/s^2 Erdbeschleunigung
ny = 1.31e-6  # m^2/s bei 10°C von Wasser


########################################################################################################################
########################################################################################################################
class CrossSection:
    """main class
    A Class that should help to generate custom cross section shapes for the SWMM software.

    Attributes:
        accuracy (int): number of decimal points to use for the export
        shape (list): descriptions of the cross section as commands in a list
        _shape_description (list): points and functions to describe the cross section
        _df_abs (pandas.DataFrame): maximum 100 points to describe the cross section in absolute values
        working_directory (str): directory where the files get saved
        unit (str): unit of entered values
        double (bool): if the cross section two separate cross sections
    """

    def __init__(self, label, description=None, height=None, width=None, working_directory='', unit=None):
        """Initialise the cross section class

        Args:
            label (str): main name/label/number of the cross section
            description (Optional[str]): optional longer name of the cross section
            height (float): absolute height of the CS
            width (Optional[float]): absolute width of the CS (optional) can be calculated
            working_directory (str): directory where the files get saved
            unit (Optional[str]): enter unit to add the unit in the plots
        """
        self.label = label
        self.description = ''
        if description is not None:
            self.description = str(description).strip()

        self._height = height
        self._width = width
        self.shape = list()
        self._shape_description = None  # NEW
        self.accuracy = 4
        self.working_directory = working_directory
        self.unit = unit
        self.double = False

        print('_' * 30)

        print(self)

        # Profile data
        self._df_abs = None

        # calculate stationary flow
        self._area_v = None
        self._r_hyd_v = None
        self._l_u_v = None
        self._r_hyd_v = None
        self._Q_v = None
        self._Q_v_params = dict(slope=None, k=None)
        self._v_v = None
        self._v_v_params = dict(slope=None, k=None)

    def __repr__(self):
        return str(self)

    def __str__(self):
        return '{}:  {}'.format(self.label, self.description)

    @property
    def height(self):
        return self._height

    @property
    def width(self):
        return self._width

    @property
    def out_filename(self):
        """
        filename of the figure/text-file to be created

        Returns:
            str: filename

        """
        return path.join(self.working_directory, str(self.label))

    def add(self, x_or_expr, y=None):
        """
        add part of cross section

        can be a:

        - function/expression
        - point (x,y) coordinates
        - boundary condition (x or y) of a surrounding function = only x or y is given and the other is :obj:`None`
        - slope (x=slope, y=unit of slope)


        Args:
            x_or_expr (Optional[float , sympy.Expr , None, CustomExpr]):

                - :obj:`float` : x coordinate or x-axis boundary or slope if any str keyword is used in argument ``y``
                - :obj:`CustomExpr` : Expression/function for the cross section part
                - :obj:`None` : if a y-axis boundary is given

            y (Optional[float,str]): y coordinate of unit of slope

                - :obj:`float` : x coordinate or x-axis boundary
                - :obj:`None` : if a x-axis boundary is given or an expression in ``x_or_expr``
                - :obj:`str` : argument x is a slope

                    - ``slope`` : ready to use slope 1 / :math:`\\Delta` y
                    - ``°slope`` : slope as an angle in degree (°)
                    - ``%slope`` : slope in percentage (%)
        """
        self._df_abs = None
        self._shape_description = None

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

    @property
    def df_abs(self, max_number_points=100):
        """create absolute point coordinates and write it into :py:attr:`~df_abs`

        To create a :obj:`pandas.DataFrame` of all the points to describe the cross section.
        This function replaces the Expressions given in :py:attr:`~add` to points with x and y coordinates
        and writes them into the :py:attr:`~df_abs` attribute.

        Args:
            max_number_points (int): number of points to describe the shape of the cross section
                                     100 is the limit of points which can be used as a SWMM shape
        """
        if self._df_abs is None:
            # number of expressions used in shape
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
                step = sum(function_steps.values()) / (max_number_points - num_points)

                min_step = 1 * 10 ** (-self.accuracy) * self.height
                if step < min_step:
                    step = min_step

            x = list()
            y = list()

            for start, end, f in self.shape_description:
                if isinstance(f, Circle):
                    # this_step = (end - start) / np.floor((end - start) / step)
                    # print(step, ' vs ', this_step)
                    nx = np.arange(start, end + step, step).clip(max=end)
                    ny = f.solve(nx)
                    x += list(nx)
                    y += list(ny)
                elif isinstance(f, Horizontal):
                    x0, y0 = f.start_point()
                    x1, y1 = f.end_point()
                    x += [x0, x1]
                    y += [y0, y1]
                else:
                    nx = np.array([start, end])
                    x += list(nx)
                    y += list(f.solve(nx))

            # the absolute points of the final shape
            df = pd.DataFrame()
            df['x'] = x
            df['y'] = y

            df = df.drop_duplicates().reset_index(drop=True)

            self._df_abs = df.astype(float).copy()

        return self._df_abs

    def get_width(self):
        """
        get absolute width of cross section

        Returns:
            float: width of cross section
        """
        if self.df_abs.empty:
            return None
        else:
            return self.df_abs['y'].max() * 2

    def set_double_cross_section(self):
        """
        make the cross section as a double section (=Doppelprofil)
        """
        self.double = True

    def check_point_cloud(self):
        """
        remove errors from the point cloud, ie.:

        - remove duplicates,
        - (if specified) remove points which overlap the overall cross section width and
        - other errors...
        """
        df = self.df_rel
        df = df.round(self.accuracy)
        df = df.drop_duplicates()

        if self.width is not None and any(df['y'] > self.width / 2):
            df['y'] = df['y'].clip_upper(self.width / 2)
            warnings.warn('had to clip the width')

        # print((np.arctan(df['x'].diff() / df['y'].diff())/ np.pi * 180).head(10))
        df = df.dropna()

        if self.double:
            df['y'] *= 2

        # delete errors
        df = df[df['x'].expanding().max() == df['x']].copy()

        # change x duplicates
        dupls = df['x'].duplicated(keep=False)
        if dupls.any():
            nx = df['x'][dupls]

            def raise_values(s):
                return s + pd.Series(index=s.index, data=range(len(s.index))) * 10 ** (-self.accuracy)

            nx = nx.groupby(nx).apply(raise_values)
            df.loc[nx.index, 'x'] = nx

        self._df_abs = (df * self.height).copy()

    @property
    def df_rel(self):
        """relative point coordinates

        convert the absolute values in the point coordinates to values relative to the cross section height

        Returns:
            pandas.DataFrame: point coordinate values relative to the cross section height
        """
        return (self.df_abs / self.height).copy()

    def make(self, show=False, plot=True):
        """
        :py:attr:`~generator` + :py:attr:`~profile_abs_plot` + :py:attr:`~input_file`

        macro function

        Args:
            show (bool):  see :py:attr:`~generator` arguments
            plot (bool): if  :py:attr:`~profile_abs_plot` should be executed
        """
        if show:
            print(self.df_abs)
        self.check_point_cloud()
        if plot:
            self.profile_abs_plot(show, file_format='pdf')
        self.input_file()

    def add_and_show(self, *args, **kwargs):
        """
        :py:attr:`~add` + :py:attr:`~generator` + :py:attr:`~profile_abs_figure`

        and print the raw shape (description of the cross section)

        macro function for jupyter example

        Args:
            *args: see :py:attr:`~add` arguments
            **kwargs: see :py:attr:`~add` keyword arguments
        """
        self.add(*args, **kwargs)
        # print('-' * 5, *self.shape, '-' * 5, sep='\n')
        print('-' * 5, *self.shape_description, '-' * 5, sep='\n')
        self.check_point_cloud()
        self.profile_abs_figure()

    def profile_rel_plot(self, auto_open=False, file_format='png'):
        """
        create a plot graphic into the :py:attr:`~working_directory` with relative dimensions

        Args:
            auto_open (bool): whether the plot should be opened after its creation
            file_format (str): file format, ie: ``png``, ``pdf``, ... (see :py:meth:`matplotlib.figure.Figure.savefig`)
        """
        ax = self.df_rel.plot(x='y', y='x', legend=False)
        ax.set_aspect('equal', 'box')
        ax.set_ylabel('rel H')
        ax.set_xlabel('B/H')
        ax.set_title('{}: {}'.format(self.label, self.description))
        fig = ax.get_figure()

        # ---------
        filename = self.out_filename + '_rel.' + file_format

        # ---------
        fig.savefig(filename)
        fig.clf()
        plt.close(fig)
        if auto_open:
            open_file(filename)

    def profile_abs_plot(self, auto_open=False, file_format='png'):
        """
        create a plot graphic into the :py:attr:`~working_directory` with absolute dimensions.

        This function uses the :py:attr:`~profile_abs_figure` -function to get a figure and saves the figure in a file.

        Args:
            auto_open (bool): whether the plot should be opened after its creation
            file_format (str): file format, ie: ``png``, ``pdf``, ... (see :py:meth:`matplotlib.figure.Figure.savefig`)
        """
        fig = self.profile_abs_figure()

        # ---------
        filename = self.out_filename + '_abs.' + file_format

        # ---------
        fig.savefig(filename)
        fig.clf()
        plt.close(fig)
        if auto_open:
            open_file(filename)

    def profile_abs_figure(self):
        """
        create a plot of the absolute dimensions

        Returns:
            matplotlib.figure.Figure: plot of the absolute dimensions
        """
        df = self.df_abs.reset_index()

        w = int(df['y'].max())
        h = int(df['x'].max())

        def custom_round(x_, base):
            return int(base * ceil(float(x_) / base))

        base = int(10 ** floor(log10(w)))
        half_base = int(base / 2)
        xlim = custom_round(w, base)
        ylim = custom_round(h, base)

        other_side = df.copy().sort_index(ascending=False)
        other_side['y'] *= -1
        df = df.append(other_side).reset_index(drop=True)
        ax = df.plot(x='y', y='x', legend=False, zorder=1000000, clip_on=False)
        ax.set_aspect('equal', 'box')
        ax.set_xticks(list(range(-xlim, xlim, base)), minor=False)
        ax.set_xticks(list(range(-xlim, xlim, half_base)), minor=True)

        ax.set_yticks(list(range(0, ylim, base)), minor=False)
        ax.set_yticks(list(range(0, ylim, half_base)), minor=True)
        # ax.set_axis_off()
        # ax.set_frame_on(False)
        # ax.axis()
        ax.tick_params(which='both', length=0, width=0, labelbottom=False, labeltop=False, labelleft=False,
                       labelright=False, bottom=False, top=False, left=False, right=False)

        ax.set_xlim(-xlim, xlim)
        ax.set_ylim(0, ylim)
        ax.grid(True)
        # ax.grid(True, which='minor', linestyle=':', linewidth=0.5)
        ax.set_xlabel(None)
        ax.set_axisbelow(True)
        # ax.set_ylabel('rel H')
        # ax.set_xlabel('B/H')

        n = self.label
        if self.label != self.description:
            n += ': {}'.format(self.description)

        ax.set_title('{}\n{:0.0f}x{:0.0f}'.format(n, h, custom_round(w * 2, half_base)) +
                     (self.unit if self.unit is not None else ''))

        fig = ax.get_figure()
        fig.tight_layout()
        return fig

    @property
    def df_rel_fill_swmm(self):
        df = self.df_rel.copy()
        df = df.iloc[1:-1].copy()
        df['y'] *= 2
        return df

    def dat_file(self):
        """
        create the EPA-SWMM Curve Data  ``.dat`` -file, which can be imported into SWMM
        The file is save into the :py:attr:`~working_directory`.
        """
        file = open(self.out_filename + '.dat', 'w+')
        file.write('EPASWMM Curve Data\n')
        dim = 'H={} {}'.format(self.height, self.unit)
        if self.width:
            dim += ', B={}'.format(self.width)
        file.write('{} - {}: {}\n'.format(self.label, self.description, dim))
        self.df_rel_fill_swmm.to_csv(file, sep=' ', index=False, header=False,
                                     float_format='%0.{}f'.format(self.accuracy))

    def inp_string(self):
        """
        create the curve data for cross section shapes in the ``.inp`` -file (SWMM-Input) format,
        which can be pasted into the input file.

        Returns:
            str: formatted text of the data
        """
        df = self.df_rel_fill_swmm

        df['name'] = path.basename(self.out_filename)

        df['shape'] = ''
        df.loc[1, 'shape'] = 'shape'

        df = df.set_index('name')

        return df[['shape', 'x', 'y']].to_string(header=False, index=True, index_names=False,
                                                 float_format='%0.{}f'.format(self.accuracy))

    def input_file(self):
        """
        create the curve data for cross section shapes in the ``.inp`` -file (SWMM-Input) format
        and save it as a separate txt-file.

        This function uses the :py:attr:`~inp_string` -function to get the string and saves the string in a file.
        The file is save into the :py:attr:`~working_directory`.
        """
        with open(self.out_filename + '_shape.txt', 'w') as f:
            f.write(self.inp_string())

    ####################################################################################################################
    # testing new functions
    ####################################################################################################################
    @property
    def shape_description(self):
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

    def b_w_t(self, hi):
        if isinstance(hi, np.ndarray):
            w = np.array([np.NaN] * hi.size)
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
        if isinstance(hi, np.ndarray):
            l = np.array([0.] * hi.size)

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

    def area_t(self, hi):
        if isinstance(hi, np.ndarray):
            a = np.array([0.] * hi.size)

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

    def r_hyd_t(self, hi):
        return self.area_t(hi) / self.l_u_t(hi)

    @property
    def area_v(self):
        if self._area_v is None:
            self._area_v = self.area_t(self.height)
        return self._area_v

    @property
    def r_hyd_v(self):
        if self._r_hyd_v is None:
            self._r_hyd_v = self.area_v / self.l_u_v
        return self._r_hyd_v

    @property
    def l_u_v(self):
        if self._l_u_v is None:
            self._l_u_v = self.l_u_t(self.height)
        return self._l_u_v

    ####################################################################################################################
    @staticmethod
    def _velocity(slope, k, r_hyd):
        """
        calculate velocity in partially filled sewer channel

        Args:
            slope (float): ablosute slope in m/m
            k (float): Betriebliche Rauhigkeit in mm
            r_hyd (float): hydraulic radius in mm

        Returns:
            float: velocity in m/s
        """
        r_hyd /= 1000
        J = slope  # / 1000
        k = k / 1000
        return (-2 * np.log10(2.51 * ny / (4 * r_hyd * sqrt(2 * g * J)) + k / (14.84 * r_hyd)) * sqrt(
            2 * g * 4 * r_hyd * J))

    def velocity_t(self, hi, slope, k):
        """
        calculate velocity in partially filled sewer channel

        Args:
            hi (float): water level = height in mm
            slope (float): ablosute slope in m/m
            k (float): Betriebliche Rauhigkeit in mm

        Returns:
            float: velocity in m/s
        """
        return self._velocity(slope, k, self.r_hyd_t(hi))

    def velocity_v(self, slope, k):
        """
        calculate velocity in partially filled sewer channel

        Args:
            slope (float): ablosute slope in m/m
            k (float): Betriebliche Rauhigkeit in mm

        Returns:
            float: full filling velocity in m/s
        """
        new_v_v_params = dict(slope=slope, k=k)
        if self._v_v is None or self._v_v_params != new_v_v_params:
            self._v_v_params = new_v_v_params
            self._v_v = self.velocity_t(self.height, slope, k)

        return self._v_v

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
            slope (float): ablosute slope in m/m
            k (float): Betriebliche Rauhigkeit in mm

        Returns:
            float: full filling flow rate in L/s
        """
        return self.velocity_v(slope, k) * self.area_v * 1.0e-6 * 1000

    ####################################################################################################################
    def h_t(self, Q_t, slope, k):
        # hi = '?'
        # self.flow_t(hi, slope, k)
        # self.flow_v(slope, k)
        from scipy.optimize import minimize_scalar
        res = minimize_scalar(lambda hi: abs(Q_t - self.flow_t(hi, slope, k)), bounds=(0, self.height),
                              method='bounded')
        return res.x


########################################################################################################################
########################################################################################################################
class CrossSectionHolding(CrossSection):
    def __init__(self, label, add_dim=False, add_dn=None, **kwargs):
        """Initialise the cross section class

        Args:
            label (str): main name/label/number of the cross section
            add_dim (bool): if the dimensions should be added to :py:attr:`~out_filename` used for the export
            add_dn (Optional[float]): if the channel dimension should be added to :py:attr:`~out_filename`
                                    used for the export enter the diameter as float
            **kwargs (object): see :py:attr:`~__init__`

        Keyword Args:
            description (Optional[str]): optional longer name of the cross section
            height (float): absolute height of the CS
            width (Optional[float]): absolute width of the CS (optional) can be calculated
            working_directory (str): directory where the files get saved
            unit (Optional[str]): enter unit to add the unit in the plots

        """
        if isinstance(label, (float, int)):
            label = '{:0.0f}'.format(label)
        else:
            label = label

        if not label.startswith('Pr_'):
            label = 'Pr_' + label

        self.add_dim = add_dim
        self.add_dn = add_dn

        CrossSection.__init__(self, label, **kwargs)

    def __str__(self):
        s = CrossSection.__str__(self)
        if self.add_dim:
            s += '  |  {:0.0f}'.format(self.height)
            if self.width:
                s += 'x{:0.0f}'.format(self.width)

        if self.add_dn:
            s += '  |  DN{:0.0f}'.format(self.add_dn)
        return s

    @property
    def out_filename(self):
        """
        filename of the figure/text-file to be created

        Returns:
            str: filename

        """
        file = path.join(self.working_directory, '{}'.format(self.label))
        if self.add_dim:
            file += '_{:0.0f}'.format(self.height)
            if self.width:
                file += 'x{:0.0f}'.format(self.width)

        if self.add_dn:
            file += '_DN{:0.0f}'.format(self.add_dn)
        return file

    ####################################################################################################################
    @classmethod
    def standard(cls, label, description, height, width=NaN, r_channel=NaN, r_roof=NaN, r_wall=NaN, slope_bench=NaN,
                 r_round=NaN, r_wall_bottom=NaN, h_bench=NaN, pre_bench=NaN, w_channel=NaN, **kwargs):
        """
        standard cross section

        Args:
            label (str): see :py:attr:`~__init__`
            description (str): see :py:attr:`~__init__`
            height (float): see :py:attr:`~__init__`
            width (float): see :py:attr:`~__init__`

            r_channel (float): radius of the dry-weather channel (=Trockenwetter Rinne)

            w_channel (float): half width of the channel, only in combination with ``r_channel`` active
            pre_bench (float): slope of the upper end of the channel in degree, only in combination with ``r_channel`` active
            r_round (float): radius of the rounding of the edges, only in combination with ``r_channel`` active
            h_bench (float): height where the bench begins, only in combination with ``r_channel`` active
            slope_bench (float): slope of the bench (=Berme) in degree, or slope of the rainwater-floor (=Regenwetterrinne)

            r_roof (float): radius of the roof (=Decke)

            r_wall (float): radius of the sidewall (=Seitenwand), only in combination with ``r_roof`` active
            r_wall_bottom (float): radius of the bottom sidewall (=untere Seitenwand), only in combination with ``r_wall`` active

            **kwargs (object): see :py:attr:`~__init__`

        Keyword Args:
            working_directory (str): directory where the files get saved
            unit (Optional[str]): enter unit to add the unit in the plots

        Returns:
            CrossSectionHolding: standard cross section

        Examples:
            see :ref:`Examples_for_standard_profiles`


        .. figure:: images/standard.gif
            :align: center
            :alt: standard cross section
            :figclass: align-center

            Standard cross section

        +---------+---------------------+
        | english | deutsch             |
        +=========+=====================+
        | channel | Trockenwetter-Rinne |
        +---------+---------------------+
        | roof    | Firste/Decke        |
        +---------+---------------------+
        | wall    | Seitenwand          |
        +---------+---------------------+
        | bench   | Berme               |
        +---------+---------------------+

        """

        # ------------------------------------------------
        cross_section = cls(label=label, description=description, height=height,
                            width=(width if notna(width) else None), **kwargs)

        # ------------------------------------------------
        # TW-Rinne
        if notna(r_channel):
            cross_section.add(Circle(r_channel, x_m=r_channel))

            # ------------------------------------------------
            if notna(pre_bench):
                cross_section.add(channel_end(r_channel, pre_bench))

                if notna(h_bench) or isna(slope_bench):
                    cross_section.add(pre_bench, '°slope')

                if notna(h_bench):
                    cross_section.add(h_bench)

            elif notna(w_channel):
                cross_section.add(None, w_channel)

            else:
                if notna(h_bench):
                    cross_section.add(h_bench)
                else:
                    cross_section.add(r_channel)
                    if isna(r_round):
                        r_round = 0
                    cross_section.add(r_channel + r_round, r_channel)

        # ------------------------------------------------
        if notna(slope_bench):
            # Berme winkel in °
            cross_section.add(slope_bench, '°slope')

        # ------------------------------------------------
        if isna(r_channel) and isna(slope_bench):
            cross_section.add(0, width / 2)

        # ------------------------------------------------
        if isna(r_roof):
            # eckige Decke
            cross_section.add(None, width / 2)
            cross_section.add(height, width / 2)

        else:
            if isna(r_wall):
                cross_section.add(None, width / 2)
                cross_section.add(height - r_roof, width / 2)
            else:
                # ------------------------------------------------
                h1 = sqrt((r_wall - r_roof) ** 2 - (r_wall - width / 2) ** 2)
                # h_middle = round(height - r_roof - h1, 8)
                h_middle = height - r_roof - h1

                # ------------------------------------------------
                if isna(r_wall_bottom):
                    cross_section.add(None, width / 2)
                    cross_section.add(h_middle, width / 2)

                else:
                    cross_section.add(Circle(r_wall_bottom, x_m=h_middle, y_m=width / 2 - r_wall_bottom))
                    cross_section.add(h_middle)

                # ------------------------------------------------
                cross_section.add(Circle(r_wall, x_m=h_middle, y_m=width / 2 - r_wall))
                cross_section.add(h_middle + h1 / (r_wall - r_roof) * r_wall)

            # ------------------------------------------------
            cross_section.add(Circle(r_roof, x_m=height - r_roof))

        # ------------------------------------------------
        return cross_section

    ####################################################################################################################
    @classmethod
    def box(cls, label, height, width, channel=None, bench=None, roof=None, rounding=0.0, **kwargs):
        """
        pre defined box (=Kasten) cross section

        see :ref:`Examples_for_box_shaped_profiles`

        Args:
            label (str): see :py:attr:`~__init__`
            height (float): see :py:attr:`~__init__`
            width (float): see :py:attr:`~__init__`
            channel (Optional[float]): diameter of the dry weather channel
            bench (Optional[float]): bench (=Berme)

                - ``''``: flache Berme
                - ``'R'``: V-förmiges Profil
                - ``'H'``: Schräge Verschneidung

            roof (Optional[float]): roof (=Decke)

                - ``''``: gerade
                - ``'B'``: Bogen
                - ``'K'``: Kreis

            rounding (Optional[float]): rounding of the edges
            **kwargs (object): see :py:attr:`~__init__`

        Keyword Args:
            description (Optional[str]): optional longer name of the cross section
            working_directory (str): directory where the files get saved
            unit (Optional[str]): enter unit to add the unit in the plots

        Returns:
            CrossSectionHolding: pre defined box (=Kasten) cross section
        """
        name = 'K'
        cross_section = cls(label=label, width=width, height=height, **kwargs)

        if isna(channel):
            channel = None
        else:
            channel = float(channel)

        if isna(bench):
            bench = ''

        if isna(roof):
            roof = ''

        if channel or bench:
            name += '.'
            bench = str(bench).strip()

            if isinstance(channel, float):
                name += '{:0.0f}'.format(channel)
                # diameter to radius
                channel /= 2
                cross_section.add(Circle(channel, x_m=channel))

            if isinstance(bench, str):
                if bench != '45':
                    name += str(bench)

                if bench == 'R':
                    cross_section.add(30, '%slope')
                    cross_section.add(None, width / 2)

                elif bench == 'H':
                    cross_section.add(channel_end(channel, 45))
                    cross_section.add(45, '°slope')
                    cross_section.add(None, width / 2)

                elif bench == '45':
                    cross_section.add(channel_end(channel, 45))
                    cross_section.add(45, '°slope')
                    cross_section.add(channel + rounding)
                    cross_section.add(channel + rounding, width / 2)

                else:
                    # Berme
                    # cross_section.add(channel, width / 2)
                    # cross_section.add(channel + rounding, width / 2)

                    if 1:
                        cross_section.add(channel_end(channel, 45))
                        cross_section.add(45, '°slope')
                    else:
                        cross_section.add(channel, channel)
                    cross_section.add(channel + rounding, None)
                    cross_section.add(5, '°slope')
                    cross_section.add(None, width / 2)

        else:
            # ebene Sohle
            cross_section.add(0, width / 2)

        if roof:
            name += '_' + str(roof)

            if roof == '':
                # gerade Decke
                cross_section.add(height, width / 2)

            elif roof == 'B':
                # Bogen-Decke
                cross_section.add(height - width * (1 - cos(radians(30))), width / 2)
                cross_section.add(Circle(width, x_m=height - width))

            elif roof == 'K':
                # Kreis Decke
                cross_section.add(height - width / 2, width / 2)
                cross_section.add(Circle(width / 2, x_m=height - width / 2))
        else:
            # gerade Decke
            cross_section.add(height, width / 2)

        if cross_section.label is None or cross_section.label == '':
            cross_section.label = name
        cross_section.check_point_cloud()
        return cross_section

    ####################################################################################################################
    @classmethod
    def box_from_string(cls, label, **kwargs):
        """
        create pre defined box (=Kasten) cross section with the string label.
        This function takes the information from the label and pass them to the :py:attr:`~box` - function.

        see :ref:`Examples_for_box_shaped_profiles`

        Args:
            label (str): see the :ref:`Examples_for_box_shaped_profiles`
            **kwargs (object): see :py:attr:`~__init__`

        Keyword Args:
            height (float): see :py:attr:`~__init__`
            width (float): see :py:attr:`~__init__`
            rounding (Optional[float]): rounding of the edges
            description (Optional[str]): optional longer name of the cross section
            working_directory (str): directory where the files get saved
            unit (Optional[str]): enter unit to add the unit in the plots

        Returns:
            CrossSectionHolding: pre defined box (=Kasten) cross section

        Examples:
            .. figure:: images/Kasten-Profile.gif
                :align: center
                :alt: Kasten-Profile
                :figclass: align-center

                Kasten-Profile
        """

        import re
        infos = re.findall(r'(K)(\.?)(\d*)([RH]?)_?([BK]?)', label)  # _(\d+)x(\d+)
        if len(infos) == 1:
            infos = infos[0]
            _, _, channel, bench, roof = infos

            if channel != '':
                channel = float(channel)
            else:
                channel = None

            bench, roof = [x_ if x_ != '' else None for x_ in (bench, roof)]

            cross_section = cls.box(label, channel=channel, bench=bench, roof=roof, **kwargs)
            return cross_section

        # --------------------------------------
        else:
            raise NotImplementedError('"{}" unknown !'.format(label))

    @classmethod
    def from_point_cloud(cls, relative_coordinates, *args, **kwargs):
        """
        get the cross sections from a point cloud where every point is relative to the lowers point in the cross section

        Args:
            relative_coordinates (pandas.Series): height-variable as index and width as values
                                                  with the origin in the lowest point of the cross section
            *args: arguments, see :py:attr:`~__init__`
            **kwargs: keyword arguments, see :py:attr:`~__init__`

        Keyword Args:
            label (str): main name/label/number of the cross section
            description (Optional[str]): optional longer name of the cross section
            height (float): absolute height of the CS
            width (Optional[float]): absolute width of the CS (optional) can be calculated
            working_directory (str): directory where the files get saved
            unit (Optional[str]): enter unit to add the unit in the plots
            add_dim (bool): if the dimensions should be added to :py:attr:`~out_filename` used for the export
            add_dn (Optional[float]): if the channel dimension should be added to :py:attr:`~out_filename`
                                    used for the export enter the diameter as float

        Returns:
            CrossSectionHolding: of the point cloud

        .. figure:: images/point_cloud.gif
            :align: center
            :alt: point cloud
            :figclass: align-center

            Point cloud
        """
        # height of the profile = maximum Y coordinate
        height_pr = relative_coordinates.index.max()

        # horizontal distance to lowest point (dry weather channel)
        # split to left and right part of the profile
        y_df = pd.concat([relative_coordinates.loc[:height_pr].rename('right'),
                          relative_coordinates.loc[height_pr:].rename('left')], axis=1)
        y_df.loc[0] = 0
        # interpolate to an even number of points on the left and right part of the profile
        y_df_filled = y_df.interpolate()

        # for debugging
        # plot of the point cloud
        # y_df_filled.plot()

        cross_section = cls(*args, **kwargs)

        if 0:
            X = 'x'  # horizontal distance to lowest point (dry weather channel)
            Y = 'y'  # vertical distance to lowest point (dry weather channel)

            df = pd.DataFrame(
                {
                    X: y_df.index,
                    Y: (y_df_filled['right'] - y_df_filled['left']) / 2
                }).reset_index(drop=True)

            cross_section._df_abs = df.copy()
            cross_section.double = False
            cross_section.check_point_cloud()

        else:
            xi = y_df.index.tolist()
            yi = ((y_df_filled['right'] - y_df_filled['left']) / 2).round(10).tolist()

            for x, y in zip(xi, yi):
                if y == 0 and x in (0, height_pr):
                    continue
                cross_section.add(x, y)

        return cross_section

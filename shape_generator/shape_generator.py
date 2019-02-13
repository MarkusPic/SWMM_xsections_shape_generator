import warnings
from math import radians, cos, ceil, log10, floor, sqrt
from numbers import Rational
from os import path
from webbrowser import open as open_file

import matplotlib.pyplot as plt
import numpy as np
import pandas
from numpy import NaN
from pandas import isna, notna
from sympy import Expr, solve, diff, Float  # , tan, cos

from .helpers import deg2slope, channel_end, circle, linear, x


########################################################################################################################
########################################################################################################################
class CrossSection(object):
    """main class
    A Class that should help to generate custom cross section shapes for the SWMM software.

    Attributes:
        accuracy (int): number of decimal points to use for the export
        shape (list): descriptions of the cross section as commands in a list
        shape_corrected (list): points and functions to describe the cross section
        df_abs (pandas.DataFrame): maximum 100 points to describe the cross section in absolute values
        working_directory (str): directory where the files get saved
        unit (str): unit of entered values
        double (bool): if the cross section two separate cross sections
    """

    def __init__(self, label, long_label=None, height=None, width=None, add_dim=False, add_dn=None,
                 working_directory='', unit=None):
        """Initialise the cross section class

        Args:
            label (str): main name/label/number of the cross section
            long_label (Optional[str]): optional longer name of the cross section
            height (float): absolute height of the CS
            width (Optional[float]): absolute width of the CS (optional) can be calculated
            add_dim (bool): if the dimensions should be added to :py:attr:`~out_filename` used for the export
            add_dn (Optional[float]): if the channel dimension should be added to :py:attr:`~out_filename`
                                    used for the export enter the diameter as float
            working_directory (str): directory where the files get saved
            unit (Optional[str]): enter unit to add the unit in the plots
        """
        if isinstance(label, (float, int)):
            self.label = '{:0.0f}'.format(label)
        else:
            self.label = label

        if not self.label.startswith('Pr_'):
            self.label = 'Pr_' + self.label

        self._name = ''
        if long_label is not None:
            self.name = str(long_label).strip()
        self.height = height
        self.width = width
        self.shape = list()
        self.shape_corrected = list()
        self.add_dim = add_dim
        self.add_dn = add_dn
        self.accuracy = 4
        self.working_directory = working_directory
        self.unit = unit
        self.double = False

        # Profile data
        self.df_abs = pandas.DataFrame()

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, value):
        self._name = value
        print('_' * 30)
        print(self.label, ' -> ', self._name)

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

    def add(self, x_or_expr, y=None):
        """
        add part of cross section

        can be a:

        - function/expression
        - point (x,y) coordinates
        - boundary condition (x or y) of a surrounding function = only x or y is given and the other is :obj:`None`
        - slope (x=slope, y=unit of slope)


        Args:
            x_or_expr (Optional[float , Expr , None , tuple]):

                - :obj:`float` : x coordinate or x-axis boundary or slope if any str keyword is used in argument ``y``
                - :obj:`Expr` : Expression/function for the cross section part
                - :obj:`None` : if a y-axis boundary is given
                - :obj:`tuple` : will be seen as the full input for this function

            y (Optional[float,str]): y coordinate of unit of slope

                - :obj:`float` : x coordinate or x-axis boundary
                - :obj:`None` : if a x-axis boundary is given or an expression in ``x_or_expr``
                - :obj:`str` : argument x is a slope

                    - ``slope`` : ready to use slope 1 / :math:`\\Delta` y
                    - ``°slope`` : slope as an angle in degree (°)
                    - ``%slope`` : slope in percentage (%)
        """
        if isinstance(x_or_expr, tuple):
            self.add(*x_or_expr)
            return
        elif isinstance(x_or_expr, Expr):
            self.shape.append(x_or_expr)
        else:
            x = x_or_expr
            if isinstance(x, (float, Rational)):
                x = float(x)
            if isinstance(y, (float, Rational)):
                y = float(y)

            if y == 'slope':
                pass
            elif y == '°slope':
                x = deg2slope(x)
                y = 'slope'
            elif y == '%slope':
                x = x / 100
                y = 'slope'

            self.shape.append((x, y))

    def check_for_slopes(self, debug=False):
        """convert slopes to points

        Use this function after adding all the necessary descriptions of the cross section with :py:attr:`~add`.
        This function converts slopes into point coordinates and specify boundary condition to (x,y) coordinates
        for Expressions and slopes.

        Args:
            debug (bool): to print debug messages during the runtime
        """
        self.shape_corrected = self.shape.copy()

        for i in range(len(self.shape_corrected)):
            if debug:
                print(i, ': ', self.shape_corrected[i], end=' ')

            # ----------------------------------------------------------------------------------------------------------
            # eine Steigung ist gegeben und wird mit End X und Y ersetzt
            if isinstance(self.shape_corrected[i], tuple) and self.shape_corrected[i][1] == 'slope':
                slope = self.shape_corrected[i][0]

                # p0: ist der letzte bekannte Punkt
                # der erste schritt ist 0,0
                if i == 0:
                    p0 = (0, 0)
                elif i == 1 and isinstance(self.shape_corrected[i - 1], tuple):
                    p0 = self.shape_corrected[i - 1]
                else:
                    p_2 = self.shape_corrected[i - 2]
                    p0_ = self.shape_corrected[i - 1]
                    x0 = p0_[0]

                    if isinstance(p_2, Expr) and isinstance(x0, float):
                        y0 = p_2.subs(x, x0)
                    elif isinstance(p_2, tuple) and isinstance(p_2[1], float):
                        y0 = p_2[1]
                    elif isinstance(p0_, tuple) and all(isinstance(i, float) for i in p0_):
                        x0, y0 = self.shape_corrected[i - 1]
                    else:
                        raise NotImplementedError

                    p0 = (x0, y0)

                # eine horizontale line
                if slope == 0:
                    pass  # do nothing
                    # self.shape_corrected[i] = p0
                    # print('--> ', fi, end=' ')
                    # print()
                    # continue

                # fi: lineare Funktion durch "p0" mit der Steigung "slope"
                fi = linear(slope, p0)
                if debug:
                    print('--> ', fi, end=' ')

                # wenn es der letzte Punkt ist, ist der nachfolgende Punkt der Scheitel
                # f2: nachfolgener Punkt oder nachfolgene Funktion
                if i == len(self.shape_corrected) - 1:
                    f2 = (float(self.height), 0.)
                else:
                    f2 = self.shape_corrected[i + 1]

                # f2: nachfolgene Funktion
                if isinstance(f2, Expr):
                    # print(type(fi))
                    if isinstance(fi, float):
                        x2 = fi
                    else:
                        try:
                            # print(f2, fi)
                            # print(solve(f2 - fi, x))
                            x2 = solve(f2 - fi, x)[0]
                        except IndexError:
                            # print(solve(diff(f2) - slope, x))
                            x2 = solve(diff(f2) - slope, x)[0]
                            if debug:
                                print('--> same slope', end=' ')

                    self.shape_corrected[i] = (float(x2), None)

                # f2: nachfolgener Punkt oder Steigung
                elif isinstance(f2, tuple):
                    # print('\n', '-'*100)
                    # print(f2)
                    # print(type(f2[0]), type(f2[1]))
                    # print()
                    if isinstance(f2[0], float):
                        if f2[1] == 'slope':
                            pass  # do nothing
                        x2 = f2[0]
                        y2 = fi.subs(x, x2)

                    elif isinstance(f2[1], float):
                        y2 = f2[1]

                        if slope != 0:
                            x2 = solve(fi - y2, x)[0]
                            self.shape_corrected[i + 1] = (float(x2), float(y2))
                        else:
                            self.shape_corrected[i + 1] = (float(p0[0]), float(y2))

                    if slope == 0:
                        self.shape_corrected[i] = (float(p0[0]), float(y2))
                    else:
                        self.shape_corrected[i] = (float(x2), float(y2))
                if debug:
                    print('--> ', self.shape_corrected[i])

            # ----------------------------------------------------------------------------------------------------------
            # Y ist gegeben und X wird ergänzt
            elif isinstance(self.shape_corrected[i], tuple) and \
                    self.shape_corrected[i][0] is None and isinstance(self.shape_corrected[i][1], float):
                if i > 0 and isinstance(self.shape_corrected[i - 1], Expr):
                    fi = self.shape_corrected[i - 1]
                    y0 = self.shape_corrected[i][1]
                    x0 = solve(fi - y0, x)[0]
                    self.shape_corrected[i] = (float(x0), y0)
                    if debug:
                        print('--> ', self.shape_corrected[i])
                else:
                    if debug:
                        print()

            else:
                if debug:
                    print()

    def create_point_cloud(self):
        """create absolute point coordinates and write it into :py:attr:`~df_abs`

        To create a :obj:`pandas.DataFrame` of all the points to describe the cross section.
        This function replaces the Expressions given in :py:attr:`~add` to points with x and y coordinates
        and writes them into the :py:attr:`~df_abs` attribute.
        """
        shape = self.shape_corrected

        # number of expressions used in shape
        num_functions = len([i for i in shape if isinstance(i, Expr)])

        # added first and last fix points
        shape = [(0., 0.)] + shape + [(self.height, 0.)]

        # if functions are used in shape
        if num_functions:
            # number of fixed points in shape
            num_points = len([i for i in shape if isinstance(i, tuple)])

            # calculate the net height of the functions.
            function_steps = {i: shape[i + 1][0] - shape[i - 1][0] for i, s in enumerate(shape) if
                              isinstance(shape[i], Expr)}
            # step size used to discretise the expressions
            step = sum(function_steps.values()) / (100 - num_points)

        # only used as an temporary variable
        # only to fill point with an expression
        is_filled = 'filled'

        # the absolute points of the final shape
        df = pandas.DataFrame(columns=['x', 'y'])

        # convert every expression to points and add it to the resulting DataFrame ``df``
        for i in range(len(shape)):
            if isinstance(shape[i], tuple):
                if isinstance(shape[i][1], type(None)):
                    continue
                if shape[i][1] == is_filled:
                    continue
                pi = shape[i]
                new = pandas.Series(list(pi), index=['x', 'y'])
                df = df.append(new, ignore_index=True)
            elif isinstance(shape[i], Expr):
                yi = shape[i]

                start = shape[i - 1][0]
                end = shape[i + 1][0]

                if start == end:
                    print('Warning: unused part of the shape detected. Ignoring this part.')
                    continue

                this_step = (end - start) / np.floor((end - start) / step)

                if isinstance(shape[i + 1][1], type(None)):
                    end += this_step
                    shape[i + 1] = (shape[i + 1][0], is_filled)

                if start == 0:
                    start += this_step
                elif not isinstance(shape[i - 1][1], type(None)):
                    start += this_step

                # x-coordinates array to discretise one expression
                xi = np.arange(start, end, this_step)

                # y-coordinates array to discretise one expression
                new = pandas.Series(xi, name='x').to_frame()
                new['y'] = np.vectorize(lambda x_i: float(yi.subs(x, Float(round(x_i, 3)))))(xi)

                df = df.append(new, ignore_index=True)

        self.df_abs = df.copy()

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
                return s + pandas.Series(index=s.index, data=range(len(s.index))) * 10 ** (-self.accuracy)

            nx = nx.groupby(nx).apply(raise_values)
            df.loc[nx.index, 'x'] = nx

        self.df_abs = (df * self.height).copy()

    @property
    def df_rel(self):
        """relative point coordinates

        convert the absolute values in the point coordinates to values relative to the cross section height

        Returns:
            pandas.DataFrame: point coordinate values relative to the cross section height
        """
        return (self.df_abs / self.height).copy()

    def generator(self, show=False):
        """
        :py:attr:`~check_for_slopes` + :py:attr:`~create_point_cloud` + :py:attr:`~check_point_cloud`

        macro function

        Args:
            show (bool): see :py:attr:`~check_for_slopes` ``debug`` - argument and print the created point cloud
        """
        self.check_for_slopes(debug=show)
        self.create_point_cloud()
        if show:
            print(self.df_abs)
        self.check_point_cloud()

    def make(self, show=False, plot=True):
        """
        :py:attr:`~generator` + :py:attr:`~profile_abs_plot` + :py:attr:`~input_file`

        macro function

        Args:
            show (bool):  see :py:attr:`~generator` arguments
            plot (bool): if  :py:attr:`~profile_abs_plot` should be executed
        """
        self.generator(show=show)

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
        print('-' * 5, *self.shape, '-' * 5, sep='\n')
        self.generator()
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
        ax.set_title('{}: {}'.format(self.label, self.name))
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
        df = self.df_abs

        w = int(df['y'].max())
        h = int(df['x'].max())

        def custom_round(x_, base):
            return int(base * ceil(float(x_) / base))

        base = int(10 ** floor(log10(w)))
        half_base = int(base / 2)
        xlim = custom_round(w, base)
        ylim = custom_round(h, base)

        other_side = df.copy().sort_values('x', ascending=False)
        other_side['y'] *= -1
        df = df.append(other_side)
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
        if self.label != 'Pr_{}'.format(self.name):
            n += ': {}'.format(self.name)

        ax.set_title('{}\n{:0.0f}x{:0.0f}'.format(n, h, custom_round(w * 2, half_base)) +
                     (self.unit if self.unit is not None else ''))
        self.cross_section_area()

        fig = ax.get_figure()
        fig.tight_layout()
        return fig

    def dat_file(self):
        """
        create the EPA-SWMM Curve Data  ``.dat`` -file, which can be imported into SWMM
        The file is save into the :py:attr:`~working_directory`.
        """
        file = open(self.out_filename + '.dat', 'w+')
        file.write('EPASWMM Curve Data\n')
        dim = 'H={} {}'.format(self.height)
        if self.width:
            dim += ', B={}'.format(self.width)
        file.write('{} - {}: {}\n'.format(self.label, self.name, dim))
        self.df_rel.iloc[1:-1].to_csv(file, sep=' ', index=False, header=False,
                                      float_format='%0.{}f'.format(self.accuracy))

    def inp_string(self):
        """
        create the curve data for cross section shapes in the ``.inp`` -file (SWMM-Input) format,
        which can be pasted into the input file.

        Returns:
            str: formatted text of the data
        """
        df = self.df_rel.copy()
        df = df.iloc[1:-1].copy()

        df['name'] = path.basename(self.out_filename)

        df['shape'] = ''
        df.loc[1, 'shape'] = 'shape'
        return df[['name', 'shape', 'x', 'y']].to_string(header=None, index=None,
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

    def cross_section_area(self):
        """
        calculate the cross section area

        Returns:
            float: area, unit depend on unit of the entered values.
        """
        # area2 = (self.df_abs['x'].diff() * self.df_abs['y'] * 2).sum()
        area = (self.df_abs['x'].diff() * (self.df_abs['y'] - self.df_abs['y'].diff() / 2) * 2).sum()
        return area

    ####################################################################################################################
    @staticmethod
    def standard(label, long_label, height, width=NaN, r_channel=NaN, r_roof=NaN, r_wall=NaN, slope_bench=NaN, r_round=NaN,
                 r_wall_bottom=NaN, h_bench=NaN, pre_bench=NaN, w_channel=NaN, add_dim=False, add_dn=False, unit=None):
        """
        standard cross section

        Args:
            label (str): see :py:attr:`~__init__`
            long_label (str): see :py:attr:`~__init__`
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

            add_dim (bool): see :py:attr:`~__init__`
            add_dn (Optional[float]): see :py:attr:`~__init__`

        Returns:
            CrossSection: standard cross section

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
        cross_section = CrossSection(label=label, long_label=long_label, height=height, width=(width if notna(width) else None),
                                     add_dim=add_dim, add_dn=add_dn, unit=unit)

        # ------------------------------------------------
        # TW-Rinne
        if notna(r_channel):
            cross_section.add(circle(r_channel, x_m=r_channel))

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
                h_middle = round(height - r_roof - h1, 8)

                # ------------------------------------------------
                if isna(r_wall_bottom):
                    cross_section.add(None, width / 2)
                    cross_section.add(h_middle, width / 2)

                else:
                    cross_section.add(circle(r_wall_bottom, x_m=h_middle, y_m=width / 2 - r_wall_bottom))
                    cross_section.add(h_middle)

                # ------------------------------------------------
                cross_section.add(circle(r_wall, x_m=h_middle, y_m=width / 2 - r_wall))
                cross_section.add(h_middle + h1 / (r_wall - r_roof) * r_wall)

            # ------------------------------------------------
            cross_section.add(circle(r_roof, x_m=height - r_roof))

        # ------------------------------------------------
        return cross_section

    ####################################################################################################################
    @staticmethod
    def box(label, height, width, channel=None, bench=None, roof=None, rounding=0.0, add_dim=True, long_label=None,
            unit=None):
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
            add_dim (bool): see :py:attr:`~__init__`
            long_label(Optional[str]): see :py:attr:`~__init__`
            unit (Optional[str]): see :py:attr:`~__init__`

        Returns:
            CrossSection: pre defined box (=Kasten) cross section
        """
        name = 'K'
        cross_section = CrossSection(label=label, long_label=long_label, width=width, height=height, add_dim=add_dim,
                                     unit=unit)

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
                cross_section.add(circle(channel, x_m=channel))

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
                cross_section.add(circle(width, x_m=height - width))

            elif roof == 'K':
                # Kreis Decke
                cross_section.add(height - width / 2, width / 2)
                cross_section.add(circle(width / 2, x_m=height - width / 2))
        else:
            # gerade Decke
            cross_section.add(height, width / 2)

        if cross_section.name is None or cross_section.name == '':
            cross_section.name = name
        cross_section.generator(show=False)
        return cross_section

    ####################################################################################################################
    @staticmethod
    def box_from_string(label, height, width, custom_label=None, unit=None):
        """
        create pre defined box (=Kasten) cross section with the string label.
        This function takes the information from the label and pass them to the :py:attr:`~box` - function.

        see :ref:`Examples_for_box_shaped_profiles`

        Args:
            label (str): see the :ref:`Examples_for_box_shaped_profiles`
            height (float): see :py:attr:`~__init__`
            width (float): see :py:attr:`~__init__`
            custom_label (str):
            unit (Optional[str]): see :py:attr:`~__init__`

        Returns:
            CrossSection: pre defined box (=Kasten) cross section

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

            cross_section = CrossSection.box(label, height=height, width=width, channel=channel, bench=bench, roof=roof,
                                             long_label=custom_label, unit=unit)
            return cross_section

        # --------------------------------------
        else:
            raise NotImplementedError('"{}" unknown !'.format(label))

    @staticmethod
    def from_point_cloud(relative_coordinates, *args, **kwargs):
        """
        get the cross sections from a point cloud where every point is relative to the lowers point in the cross section

        Args:
            relative_coordinates (pandas.Series): height-variable as index and width as values
                                                  with the origin in the lowest point of the cross section
            *args: arguments, see :py:attr:`~__init__`
            **kwargs: keyword arguments, see :py:attr:`~__init__`

        Returns:
            CrossSection: of the point cloud

        .. figure:: images/point_cloud.gif
            :align: center
            :alt: point cloud
            :figclass: align-center

            Point cloud
        """
        X = 'x'  # horizontal distance to lowest point (dry weather channel)
        Y = 'y'  # vertical distance to lowest point (dry weather channel)

        # height of the profile = maximum Y coordinate
        height_pr = relative_coordinates.index.max()

        # horizontal distance to lowest point (dry weather channel)
        # split to left and right part of the profile
        y_df = pandas.concat([relative_coordinates.loc[:height_pr].rename('right'),
                              relative_coordinates.loc[height_pr:].rename('left')], axis=1)
        y_df.loc[0] = 0
        # interpolate to an even number of points on the left and right part of the profile
        y_df_filled = y_df.interpolate()

        # for debugging
        # plot of the point cloud
        # y_df_filled.plot()

        df = pandas.DataFrame(
            {
                X: y_df.index,
                Y: (y_df_filled['right'] - y_df_filled['left']) / 2
            })

        cross_section = CrossSection(*args, **kwargs)
        cross_section.df_abs = df.copy()
        cross_section.check_point_cloud(double=False)
        return cross_section

import warnings
from sympy import symbols, Expr, sqrt, Symbol, solve, diff, Float  # , tan, cos
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from webbrowser import open as show_file
from os import path
from math import radians, cos, atan, tan, sin, acos, degrees
from numpy import NaN
from pandas import isna, notna, read_excel, read_csv
from numbers import Rational

from .helpers import deg2slope, channel_end, circle, linear, x


########################################################################################################################
########################################################################################################################
class Profile(object):
    """A Class that should help to generate custom cross section shapes for the SWMM software."""

    def __init__(self, number, name=None, height=None, width=None, add_dim=False, add_dn=False, working_directory=''):
        """

        Args:
            number (str):
            name (str):
            height (float):
            width (float):
            add_dim (bool):
            add_dn (bool):
            working_directory (str):
        """
        if isinstance(number, (float, int)):
            self.number = '{:0.0f}'.format(number)
        else:
            self.number = number

        if not self.number.startswith('Pr_'):
            self.number = 'Pr_' + self.number

        self._name = ''
        if name is not None:
            self.name = str(name).strip()
        self.height = height
        self.width = width
        self.shape = list()
        self.shape_corrected = list()
        self.add_dim = add_dim
        self.add_dn = add_dn
        self.accuracy = 4
        self.out_path = working_directory

        # Profile data
        self.df_abs = pd.DataFrame()

    @property
    def name(self):
        return self._name

    def get_width(self):
        if self.df_abs.empty:
            return None
        else:
            return self.df_abs['y'].max() * 2

    @property
    def out_filename(self, long=False):
        if long:
            file = path.join(self.out_path, '{}_{}'.format(self.number, self.name))
        else:
            file = path.join(self.out_path, '{}'.format(self.number))
        if self.add_dim:
            file += '_{:0.0f}'.format(self.height)
            if self.width:
                file += 'x{:0.0f}'.format(self.width)

        if self.add_dn:
            file += '_DN{:0.0f}'.format(self.add_dn)
        return file

    @name.setter
    def name(self, value):
        self._name = value
        print('_' * 30)
        print(self.number, ' -> ', self._name)

    def add(self, x_or_expr, y=None):
        """
        add part of cross section
        can be a:
        - function/expression
        - point (x,y)
        - boundary condition (x or y) of a surrounding function
        - slope (x=slope, y=unit of slope)

        :param x_or_expr:
        :type x_or_expr: float | Expr | None
        :param y: y coordinate of unit of slope
        :type y: float | str
        """
        if isinstance(x_or_expr, tuple):
            self.add(*x_or_expr)
            # self.shape.append((float(i) for i in x_or_expr))
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

    def check_for_slopes(self, show=False):
        """
        convert slopes into point coordinates
        and specify boundary condition to (x,y) coordinates

        :param show: for debugging
        :type show: bool
        """
        self.shape_corrected = self.shape.copy()

        for i in range(len(self.shape_corrected)):
            if show:
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
                if show:
                    print('--> ', fi, end=' ')

                # wenn es der letzte Punkt ist, ist der nachfolgende Punkt der Scheitel
                # f2: nachfolgener Punkt oder nachfolgene Funktion
                if i == len(self.shape_corrected) - 1:
                    f2 = (self.height, 0)
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
                            if show:
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
                if show:
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
                    if show:
                        print('--> ', self.shape_corrected[i])
                else:
                    if show:
                        print()

            else:
                if show:
                    print()

    def create_point_cloud(self):
        """
        create a point cloud out of the functions

        :rtype: pd.DataFrame
        :return:
        """
        # print(*self.p, sep='\n')
        try:
            # print('*' * 15)
            # print(*self.p, sep='\n')
            (shape, width, height, name, number, add_dim) = (
                self.shape_corrected, self.width, self.height, self.name, self.number, self.add_dim)

            num_funktions = len([i for i in shape if isinstance(i, Expr)])

            # if height is None:
            #     height = 10000

            shape = [(0, 0)] + shape + [(height, 0)]

            if num_funktions:
                num_points = len([i for i in shape if isinstance(i, tuple)])
                # num_shapes = len(shape)
                # if num_shapes == 1:
                funktion_steps = {i: shape[i + 1][0] - shape[i - 1][0] for i, s in enumerate(shape) if
                                  isinstance(shape[i], Expr)}
                step = sum(funktion_steps.values()) / (100 - num_points)

            is_filled = 'filled'

            df = pd.DataFrame(columns=['x', 'y'])
            for i in range(len(shape)):
                if isinstance(shape[i], tuple):
                    if isinstance(shape[i][1], type(None)):
                        continue
                    if shape[i][1] == is_filled:
                        continue
                    pi = shape[i]
                    new = pd.Series(list(pi), index=['x', 'y'])
                    df = df.append(new, ignore_index=True)
                elif isinstance(shape[i], Expr):
                    yi = shape[i]

                    start = shape[i - 1][0]
                    end = shape[i + 1][0]

                    this_step = (end - start) / np.floor((end - start) / step)

                    if isinstance(shape[i + 1][1], type(None)):
                        end += this_step
                        shape[i + 1] = (shape[i + 1][0], is_filled)

                    if start == 0:
                        start += this_step
                    elif not isinstance(shape[i - 1][1], type(None)):
                        start += this_step

                    try:
                        xi = np.arange(start, end, this_step)
                    except ValueError:
                        print(i, yi)
                        print(start, end, this_step, step)
                        raise ValueError

                    if not xi.size:
                        print(i, yi)
                        print(start, end, this_step, step)

                    new = pd.Series(xi, name='x').to_frame()
                    try:
                        new['y'] = np.vectorize(lambda x_i: float(yi.subs(x, Float(round(x_i, 3)))))(xi)
                    except ValueError:
                        print()
                        print(yi)
                        print()
                        print(xi)
                        exit()

                    df = df.append(new, ignore_index=True)

        except:
            print()
            print('#' * 30)
            print(*self.shape, sep='\n')
            print('#' * 30)
            raise ArithmeticError

        return df

    def check_point_cloud(self, df, double=False):
        """
        remove errors from point cloud and create the

        :param df:
        :param double: für Doppelfrofile
        :return:
        """
        self.df_abs = df.astype(float).copy()
        df = self.df_rel
        df = df.round(self.accuracy)
        df = df.drop_duplicates()

        if self.width is not None and any(df['y'] > self.width / 2):
            df['y'] = df['y'].clip_upper(self.width / 2)
            warnings.warn('had to clip the width')

        # print((np.arctan(df['x'].diff() / df['y'].diff())/ np.pi * 180).head(10))
        df = df.dropna()

        if double:
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

        self.df_abs = (df * self.height).copy()

    @property
    def df_rel(self):
        return (self.df_abs / self.height).copy()

    def generator(self, double=False, show=False):
        self.check_for_slopes(show=show)
        df = self.create_point_cloud()
        if show:
            print(df)
        self.check_point_cloud(df, double)

    def make(self, double=False, show=False, plot=True):
        """
        just a macro

        :type plot: bool
        :param bool double:
        :param bool show:
        """
        self.generator(double=double, show=show)

        if plot:
            self.profile_abs_plot(show, file_format='pdf')
        # self.dat_file()
        self.input_file()

    def add_and_show(self, *args, **kwargs):
        #  for jupyter
        self.add(*args, **kwargs)
        print(self.shape)
        self.generator(show=False)
        self.profile_abs_figure()

    def profile_rel_plot(self, show=False, file_format='png'):
        """
        create a png plot into

        :type file_format: str
        :param show: open the plot file in a viewer
        """
        ax = self.df_rel.plot(x='y', y='x', legend=False)
        ax.set_aspect('equal', 'box')
        ax.set_ylabel('rel H')
        ax.set_xlabel('B/H')
        ax.set_title('{}: {}'.format(self.number, self.name))

        if self.height < 10:
            self.height = int(self.height * 1000)

        filename = self.out_filename + '.' + file_format

        fig = ax.get_figure()
        fig.savefig(filename)
        # print(filename)
        fig.clf()
        plt.close(fig)
        if show:
            show_file(filename)

    def profile_abs_plot(self, show=False, file_format='pdf'):
        fig = self.profile_abs_figure()
        filename = self.out_filename + '.' + file_format

        fig.savefig(filename)
        # print(filename)
        fig.clf()
        plt.close(fig)
        if show:
            show_file(filename)

    def profile_abs_figure(self):
        """
        create a png plot into

        :type file_format: str
        :param show: open the plot file in a viewer
        """

        if self.number == 'Pr_18':
            df = self.df_rel * 1950
        elif self.number == 'Pr_66':
            df = self.df_rel * 2100
        else:
            df = self.df_abs

        w = int(df['y'].max())
        h = int(df['x'].max())

        from math import ceil
        def custom_round(x, base):
            res = int(base * ceil(float(x) / base))
            return res

        xlim = custom_round(w, 100)
        ylim = custom_round(h, 100)

        other_side = df.copy().sort_values('x', ascending=False)
        other_side['y'] *= -1
        df = df.append(other_side)
        ax = df.plot(x='y', y='x', legend=False, zorder=1000000, clip_on=False)
        ax.set_aspect('equal', 'box')
        ax.set_xticks(list(range(-xlim, xlim, 100)), minor=False)
        ax.set_xticks(list(range(-xlim, xlim, 50)), minor=True)

        ax.set_yticks(list(range(0, ylim, 100)), minor=False)
        ax.set_yticks(list(range(0, ylim, 50)), minor=True)
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

        n = self.number
        if self.number != 'Pr_{}'.format(self.name):
            n += ': {}'.format(self.name)

        ax.set_title('{}\n{:0.0f}x{:0.0f}mm'.format(n, h, custom_round(w * 2, 50)))
        self.cross_section_area()

        if self.height < 10:
            self.height = int(self.height * 1000)
        fig = ax.get_figure()
        fig.tight_layout()
        return fig

    def dat_file(self):
        """
        create the EPASWMM Curve Data file
        to import the file in SWMM
        """
        file = self.out_filename

        csv = open(file + '.dat', 'w+')
        csv.write('EPASWMM Curve Data\n')
        dim = 'H={}'.format(self.height)
        if self.width:
            dim += ', B={}'.format(self.width)
        csv.write('{} - {}: {}\n'.format(self.number, self.name, dim))
        self.df_rel.iloc[1:-1].to_csv(csv, sep=' ', index=False, header=False,
                                      float_format='%0.{}f'.format(self.accuracy))

    def inp_string(self):
        df = self.df_rel.copy()
        df = df.iloc[1:-1].copy()

        df['name'] = path.basename(self.out_filename)

        df['shape'] = ''
        df.loc[1, 'shape'] = 'shape'
        return df[['name', 'shape', 'x', 'y']].to_string(header=None, index=None,
                                                         float_format='%0.{}f'.format(self.accuracy))

    def input_file(self):
        """
        create the profile table for the ".inp"-SWMM-file as a seperate txt-file
        """
        with open(self.out_filename + '_shape.txt', 'w') as f:
            f.write(self.inp_string())

    def cross_section_area(self):
        """
        calculate the cross section a

        :rtype: float
        :return: area in m²
        """
        # area2 = (self.df_abs['x'].diff() * self.df_abs['y'] * 2).sum()
        area = (self.df_abs['x'].diff() * (self.df_abs['y'] - self.df_abs['y'].diff() / 2) * 2).sum() * 1e-6
        return area

    ####################################################################################################################
    @staticmethod
    def standard(no, name, height, width=NaN, r_channel=NaN, r_roof=NaN, r_wall=NaN, slope_bench=NaN, r_round=NaN,
                 r_wall_bottom=NaN, h_bench=NaN, pre_bench=NaN, w_channel=NaN, add_dim=False, add_dn=False):
        """
        Der Querschnitt des Standard Profils

        :type add_dn: bool
        :type add_dim: bool
        :param int | str no: number
        :param str name: label of the profile
        :param float height:
        :param float width:
        :param float r_channel: radius
        :param float r_roof: radius
        :param float r_wall: radius
        :param float slope_bench: slope in degree
        :param float r_round: radius
        :param float r_wall_bottom: radius
        :param float h_bench: height
        :param float pre_bench: fist bench in degree
        :param float w_channel: width

        :rtype: Profile
        """

        # ------------------------------------------------
        cross_section = Profile(number=no, name=name, height=height, width=(width if notna(width) else None),
                                add_dim=add_dim, add_dn=add_dn)

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
    def box(no, height, width, channel=None, bench=None, roof=None, rounding=0.0, add_dim=True,
            custom_label=None):
        """
        Kasten

        :type add_dim: bool
        :type custom_label: str
        :type no: profile number
        :param float height: in [mm]
        :param float width: in [mm]
        :param float channel: diameter in [mm]
        :param str bench: ''=flache Berme <|> 'R'=V-förmiges Profil <|> 'H'=Schräge Verschneidung
        :param str roof: ''=gerade <|> 'B'= Bogen <|>  'K'=Kreis
        :type rounding: float

        :rtype: Profile
        :return:
        """
        name = 'K'
        cross_section = Profile(number=no, name=custom_label, width=width, height=height, add_dim=add_dim)

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
                name += '{:0.0f}'.format(channel / 10)
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
    def box_from_string(label, height, width, custom_label=None):
        """

        :type custom_label: str
        :type width: float
        :type height: float
        :param label: label of box profile

        :return:
        :rtype: Profile
        """

        import re
        infos = re.findall(r'(K)(\.?)(\d*)([RH]?)_?([BK]?)', label)  # _(\d+)x(\d+)
        if len(infos) == 1:
            infos = infos[0]
            _, _, channel, bench, roof = infos

            if channel != '':
                channel = float(channel) * 10  # from cm to mm
            else:
                channel = None

            bench, roof = [x if x != '' else None for x in (bench, roof)]

            cross_section = Profile.box(label, height=height, width=width, channel=channel, bench=bench, roof=roof,
                                        custom_label=custom_label)
            return cross_section

        # --------------------------------------
        else:
            raise NotImplementedError('"{}" unknown !'.format(label))

    @staticmethod
    def from_point_cloud(excel_filename='/home/markus/Downloads/Haltung 6560086.xlsx',
                         number=950, name='Profil 950', height=2250, width=1800, add_dim=True, add_dn=False):
        """

        :param excel_filename:
        :param number:
        :param name:
        :param height:
        :param width:
        :param add_dim:
        :param add_dn:

        :rtype: Profile
        :return:
        """
        X = 'x'  # horizontal distance to lowest point (dry weather channel)
        Y = 'y'  # vertical distance to lowest point (dry weather channel)

        # distances in meter
        if excel_filename.endswith('.csv'):
            coordinates = \
            read_csv(excel_filename, header=0, usecols=[0, 1], names=[X, Y]).mul(1000).round(0).set_index(Y)[
                X]
        elif excel_filename.endswith('.xlsx'):
            coordinates = \
                read_excel(excel_filename, skiprows=3, header=None, usecols=[0, 1], names=[X, Y]).mul(1000).set_index(
                    Y)[X]
        else:
            raise NotImplementedError

        # height of the profile = maximum Y coordinate
        height_pr = coordinates.index.max()

        # horizontal distance to lowest point (dry weather channel)
        # split to left and right part of the profile
        y_df = pd.concat([coordinates.loc[:height_pr].rename('right'),
                          coordinates.loc[height_pr:].rename('left')], axis=1)
        y_df.loc[0] = 0
        # interpolate to an even number of points on the left and right part of the profile
        y_df_filled = y_df.interpolate()

        # for debugging
        # plot of the point cloud
        y_df_filled.plot()

        # s = (y_df_filled['right'] - y_df_filled['left']) / 2
        # df = pd.DataFrame()
        # df[X] = s.index.values
        # df[Y] = s.values

        df = pd.DataFrame(
            {
                X: y_df.index,
                Y: (y_df_filled['right'] - y_df_filled['left']) / 2
            })

        cross_section = Profile(number=number, name=name, height=height, width=width, add_dim=add_dim, add_dn=add_dn)
        cross_section.df_abs = df
        cross_section.check_point_cloud(df, double=False)
        return cross_section

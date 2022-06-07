from .helpers import Circle, Slope, PowerExpr
from .shape_generator import CrossSection


class Circular(CrossSection):
    def __init__(self, diameter, label=None):
        r = diameter/2
        height = diameter
        width = diameter

        if label is None:
            label = f'DN {diameter:0.0f}'

        CrossSection.__init__(self, label=label, width=width, height=height)
        self.add(Circle(r, x_m=r))


class CircularFilled(CrossSection):
    def __init__(self, diameter, filled_depth, label=None):
        r = diameter/2
        height = diameter
        width = diameter

        if label is None:
            label = f'DN {diameter:0.0f} ({filled_depth} sediments)'

        width_bottom = (r**2 - (r-filled_depth)**2)**(1/2)

        CrossSection.__init__(self, label=label, width=width, height=height)
        self.add(0, width_bottom)
        self.add(Circle(r, x_m=r-filled_depth))


class Egg(CrossSection):
    def __init__(self, height, label=None):
        r = height/3
        R = 3 * r
        rho = r / 2
        # height = r * 3
        width = r * 2
        # h1 = rho - (r + rho) / (R - rho) * rho
        h1 = 0.2 * r

        if label is None:
            label = f'Ei {width:0.0f}/{height:0.0f}'

        CrossSection.__init__(self, label=label, width=width, height=height)
        self.add(Circle(rho, x_m=rho))
        self.add(h1)
        self.add(Circle(R, x_m=2 * r, y_m=-(R - r)))
        self.add(2 * r)
        self.add(Circle(r, x_m=2 * r))


class RectangularClosed(CrossSection):
    def __init__(self, height, width, label=None):
        CrossSection.__init__(self, label=label, width=width, height=height)
        self.add(0, width/2)
        self.add(height, width/2)


class RectangularOpen(RectangularClosed):
    def __init__(self, height, width, label=None):
        RectangularClosed.__init__(self, label=label, width=width, height=height)


class RectangularRound(CrossSection):
    def __init__(self, height, width, radius_bottom, label=None):
        CrossSection.__init__(self, label=label, width=width, height=height)

        if radius_bottom >= width/2:
            self.add(Circle(radius_bottom, x_m=radius_bottom))
            height_round = radius_bottom - (radius_bottom ** 2 - (width/2) ** 2) ** (1 / 2)
            self.add(height_round)
            # self.add(None, width/2)
            self.add(height, width / 2)


class RectangularTriangular(CrossSection):
    def __init__(self, height, width, height_triangular, label=None):
        CrossSection.__init__(self, label=label, width=width, height=height)
        self.add(height_triangular, width / 2)
        self.add(height, width / 2)


class Triangular(CrossSection):
    def __init__(self, height, width_max, label=None):
        CrossSection.__init__(self, label=label, width=width_max, height=height)
        self.add(height, width_max / 2)


class Power(CrossSection):
    def __init__(self, height, width_max, exponent, label=None):
        CrossSection.__init__(self, label=label, width=width_max, height=height)
        self.add(PowerExpr(exponent=exponent, p=width_max/(2*height**(1/exponent))))


class Parabolic(Power):
    def __init__(self, height, width_max, label=None):
        Power.__init__(self, label=label, width_max=width_max, exponent=2, height=height)


# class Trapezoidal(CrossSection):
#     def __init__(self, height, width_base, slope_left, slope_right, label=None):
#         CrossSection.__init__(self, label=label, width=width, height=height)
#         self.add(0, width_base/2)
#         self.add(Slope())
#         self.add(height, width/2)

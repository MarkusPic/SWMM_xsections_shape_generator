from . import Circle, CrossSection


class EggSection(CrossSection):
    """
    egg shaped cross section

    .. figure:: images/ei.gif
            :align: center
            :alt: egg
            :figclass: align-center

            Egg Section (DWA-A 110, 2006)
    """
    def __init__(self, r, label=None, description=None):
        """init
        egg shaped cross section

        Args:
            r (float): radius of the egg
            label (str): name/label/number of the cross section; dafault = "Ei <width>/<height>"
            description (str): optional description of the cross section
        """
        R = 3 * r
        roh = r / 2
        height = r * 3
        width = r * 2
        # h1 = roh - (r + roh) / (R - roh) * roh
        h1 = 0.2 * r

        if label is None:
            label = 'Ei {:0.0f}/{:0.0f}'.format(width, height)

        CrossSection.__init__(self, label=label, description=description, width=width, height=height)
        self.add(Circle(roh, x_m=roh))
        self.add(h1)
        self.add(Circle(R, x_m=2 * r, y_m=-(R - r)))
        self.add(2 * r)
        self.add(Circle(r, x_m=2 * r))


class CircleSection(CrossSection):
    """
    circle cross section

    .. figure:: images/kreis1.gif
            :align: center
            :alt: circle
            :figclass: align-center

            Circle Section (DWA-A 110, 2006)
    """
    def __init__(self, r, label=None, description=None):
        """init
        circle cross section

        Args:
            r (float): radius of the circle
            label (str): name/label/number of the cross section; dafault = "DN <diameter>"
            description (str): optional description of the cross section
        """
        d = 2 * r
        height = d
        width = d

        if label is None:
            label = 'DN {:0.0f}'.format(d)

        CrossSection.__init__(self, label=label, description=description, width=width, height=height)
        self.add(Circle(r, x_m=r))

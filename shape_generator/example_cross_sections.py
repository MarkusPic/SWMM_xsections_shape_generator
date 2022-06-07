from .helpers import Circle
from .shape_generator import CrossSection
import json
import os

import numpy as np


class EggSection(CrossSection):
    """
    egg shaped cross section

    .. figure:: images/ei.gif
            :align: center
            :alt: egg
            :figclass: align-center

            Egg Section (DWA-A 110, 2006)
    """

    def __init__(self, r, label=None):
        """init
        egg shaped cross section

        Args:
            r (float): radius of the egg
            label (str): name/label/number of the cross section; dafault = "Ei <width>/<height>"
            description (str): optional description of the cross section
        """
        R = 3 * r
        rho = r / 2
        height = r * 3
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


class CircleSection(CrossSection):
    """
    circle cross section

    .. figure:: images/kreis1.gif
            :align: center
            :alt: circle
            :figclass: align-center

            Circle Section (DWA-A 110, 2006)
    """

    def __init__(self, r, label=None):
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
            label = f'DN {d:0.0f}'

        CrossSection.__init__(self, label=label, width=width, height=height)
        self.add(Circle(r, x_m=r))


# -------------------------------------------------
# Cross-sections pre-defined in SWMM
SWMM_STD_CROSS_SECTION_CURVES = None


def _load_swmm_std_cross_section_curves():
    global SWMM_STD_CROSS_SECTION_CURVES
    if SWMM_STD_CROSS_SECTION_CURVES is None:
        SWMM_STD_CROSS_SECTION_CURVES = json.load(
            open(os.path.join(os.path.dirname(__file__), 'swmm_std_cross_section_curves.json'), 'r'))


width_max_factor = {
    'EGG'         : 2 / 3,
    'GOTHIC'      : 0.84,
    'CATENARY'    : 0.9,
    'BASKETHANDLE': 0.944,
    'SEMICIRCULAR': 1.64,
}


def swmm_std_cross_sections(shape, height=1, width=None, label=None):
    """
    get a SWMM pre-defined cross-section

    Args:
        shape (str): name of the cross-section. one of:
            - CIRCULAR
            - EGG
            - HORSESHOE
            - GOTHIC
            - CATENARY
            - SEMIELLIPTICAL
            - BASKETHANDLE
            - SEMICIRCULAR

            - ARCH
            - HORIZ_ELLIPSE
            - VERT_ELLIPSE

        height (float): height of the cross-section
        width (float | Optional): width of the cross-section

    Returns:
        CrossSection:
    """
    global SWMM_STD_CROSS_SECTION_CURVES

    if SWMM_STD_CROSS_SECTION_CURVES is None:
        _load_swmm_std_cross_section_curves()

    if shape not in SWMM_STD_CROSS_SECTION_CURVES:
        return
    rel_with = np.array(SWMM_STD_CROSS_SECTION_CURVES[shape])*height*width_max_factor.get(shape, 1)
    rel_heights = np.linspace(0, 1, len(rel_with))*height

    cross_section = CrossSection(label, height=height)

    if width is not None:
        factor_width = width/height
    else:
        factor_width = 1

    for x, y in zip(rel_heights.round(8), rel_with.round(8)):
        if (y == 0) and (x in (0, 1)):
            continue
        cross_section.add(x, (y * factor_width) / 2)

    return cross_section

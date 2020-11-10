from shape_generator import CrossSection
from swmm_api.input_file.inp_sections import Curve
import numpy as np


def convert_shape_generator_to_curve(gen: CrossSection) -> Curve:
    """create a swmm curve object with the data of the swmm-shape_generator-CrossSection"""
    x, y = gen.get_points()
    # [1:-1] without first and last point / not needed in swmm
    height = np.array(x[1:-1]) / gen.height
    area = np.array(y[1:-1]) / gen.height * 2
    return Curve(Name=gen.label, Type=Curve.TYPES.SHAPE, points=[[float(h), float(a)] for h, a in zip(height, area)])


def from_swmm_shape(curve: Curve, height=100, *args, **kwargs) -> CrossSection:
    """
    create a object with the data of the swmm curve data as relative coordinates

    Args:
        curve (Curve): Curve object of the CURVES section in the inp-data file
        height (float): absolute height of the CS
        *args: arguments, see :attr:`CrossSection.__init__`
        **kwargs: keyword arguments, see :attr:`CrossSection.__init__`

    Returns:
        CrossSection: of the shape coordinates
    """
    cross_section = CrossSection(curve.Name, height=height, *args, **kwargs)
    for x, y in curve.points:
        cross_section.add(x*height, y*height / 2)
    return cross_section

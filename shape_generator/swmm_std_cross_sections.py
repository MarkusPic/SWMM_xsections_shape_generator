import json
import os.path

import numpy as np

from shape_generator import CrossSection

SWMM_STD_CROSS_SECTION_CURVES = json.load(open(os.path.join(os.path.dirname(__file__), 'swmm_std_cross_section_curves.json'), 'r'))


def swmm_std_cross_sections(shape, height=1):
    if shape not in SWMM_STD_CROSS_SECTION_CURVES:
        return
    rel_with = SWMM_STD_CROSS_SECTION_CURVES[shape]
    rel_heights = np.linspace(0, 1, len(rel_with))
    cross_section = CrossSection(shape, height=height)
    for x, y in zip(rel_heights, rel_with):
        if (y == 0) and (x in (0, 1)):
            continue
        cross_section.add(x, y / 2)

    return cross_section

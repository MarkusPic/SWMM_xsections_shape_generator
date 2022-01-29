Main Function
-----------------

CrossSection
================

.. currentmodule:: shape_generator.shape_generator

.. rubric:: Class

.. autosummary::
    CrossSection

.. rubric:: Methods

.. currentmodule:: shape_generator.shape_generator.CrossSection

.. autosummary::
    __init__
    identifier
    height
    width
    get_width
    add
    set_double_cross_section
    shape_description
    get_points

.. rubric:: Figures

.. autosummary::
    profile_axis
    profile_figure

.. rubric:: Shape parameters

.. autosummary::
    b_w_t
    l_u_t
    area_t
    r_hyd_t
    area_v
    r_hyd_v
    l_u_v
    velocity_t
    velocity_v
    flow_t
    flow_v
    h_t

.. rubric:: swmm_api Functions

.. autosummary::
    from_curve
    to_curve

CrossSectionHolding
=======================

.. rubric:: Class

.. currentmodule:: shape_generator.shape_generator_holding

.. autosummary::
    CrossSectionHolding

.. autoclass:: shape_generator.shape_generator.CrossSection
    :members:

.. currentmodule:: shape_generator.shape_generator_holding.CrossSectionHolding

.. autosummary::
    __init__

.. rubric:: Pre defined Cross Sections

.. autosummary::
    standard
    box
    box_from_string
    from_point_cloud

.. autoclass:: shape_generator.shape_generator_holding.CrossSectionHolding
    :members:

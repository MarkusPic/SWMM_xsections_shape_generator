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
    height
    width
    add
    get_points
    df_abs
    get_width
    set_double_cross_section
    check_point_cloud
    df_rel
    df_rel_fill_swmm
    shape_description

.. rubric:: Macros

.. autosummary::
    add_and_show

.. rubric:: Text-files / Figures

.. autosummary::
    profile_rel_plot
    profile_abs_plot
    profile_abs_figure
    dat_file
    inp_string
    input_file

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

CrossSectionHolding
=======================

.. rubric:: Class

.. currentmodule:: shape_generator.shape_generator

.. autosummary::
    CrossSectionHolding

.. currentmodule:: shape_generator.shape_generator.CrossSectionHolding

.. autosummary::
    __init__

.. rubric:: Pre defined Cross Sections

.. autosummary::
    standard
    box
    box_from_string
    from_point_cloud

CrossSectionMisc
=======================

.. rubric:: Class

.. currentmodule:: shape_generator.shape_generator

.. autosummary::
    CrossSectionMisc

.. currentmodule:: shape_generator.shape_generator.CrossSectionMisc

.. autosummary::
    __init__

.. rubric:: Pre defined Cross Sections

.. autosummary::
    from_swmm_shape

.. autoclass:: shape_generator.shape_generator.CrossSection
    :members:


.. autoclass:: shape_generator.shape_generator.CrossSectionHolding
    :members:


.. autoclass:: shape_generator.shape_generator.CrossSectionMisc
    :members:

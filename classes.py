from collections import namedtuple

Parameters = namedtuple('Parameters', [
    'L_w',
    'V_w',
    'rho_l',
    'rho_v',
    'sigma',
    'turn_down',
    't_shell',
    'l_calm',
    'Nplates',
])

Design = namedtuple('Design', [
    'd_col',
    'theta',
    'd_hole',
    'l_pitch',
    'l_spacing',
    'h_weir',
    'h_ap',
    't_plate',
])


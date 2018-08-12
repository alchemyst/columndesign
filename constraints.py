import numpy

def table_K1(Flv, l_spacing):
    """Figure 11-27
    
    Regression by Darren Roos
    """
    c, b, a, d = 0.6405720204, -0.2039915816, -0.5905334053, -1.8291520294
    return 10**(a*numpy.log10(Flv) + b*(numpy.log10(Flv))**2 + c*l_spacing + d)


def table_psi(Flv, flood):
    """Figure 11-29

    Regression by Darren Roos
    """
    c, b, a, d = -0.0137855596, 0.0057185289, -0.0722581746, -3.3343391415
    return 10**(a*numpy.log10(Flv) + b*flood + c*flood*numpy.log10(Flv) + d)


def table_K2(h):
    """Figure 11-30
    
    Regression by Carl Sandrock"""
    a, b, c, d = 26.52258403,  0.76197575, 13.23116619, 33.1867269
    return a + b*numpy.sqrt(numpy.abs(h - c)) - (h - c)/d


def table_C0(A_ratio, P_ratio):
    """Figure 11-34

    Regression by Darren Roos
    """
    b, a, c = 0.1357643198, 0.0079299573, 0.620738345
    return a*A_ratio + b*P_ratio**2 + c


CONSTRAINT_NAMES = ["Flooding", 
                    "Entrainment", 
                    "Weeping", 
                    "Downcomer backup", 
                    "Residence time"]

def columnconstraints(parameters, design):
    V, L, rho_v, rho_l, mu, sigma, turn_down, t_plate, l_calm = parameters
    d_col, d_hole, theta, l_spacing, l_pitch, h_weir, h_ap = design

    # Geometry calculations
    l_weir = d_col * numpy.sin(theta / 2)  # calculated using trig

    h = d_col / 2 * numpy.cos(theta / 2)
    h1 = h - l_calm
    d1 = d_col - 2 * l_calm
    theta1 = 2*numpy.arccos(h1 / (d1 / 2))
    A1 = (theta1 - numpy.sin(theta1)) / 2 * (d1 / 2) ** 2

    # Areas
    A_col = numpy.pi / 4 * d_col ** 2
    A_h = numpy.pi / 4 * d_hole ** 2
    A_down = (theta - numpy.sin(theta)) / 2 * (d_col / 2) ** 2
    # A_active = A_col - 2 * A_down
    A_p = numpy.pi / 4 * d1 ** 2 - 2 * A1
    A_hs = (A_h / 2) / (numpy.sqrt(3) / 4 * l_pitch ** 2) * A_p
    A_net = A_col - A_down

    def _constraints(V_w, L_w):
        # Flooding
        u = V_w / rho_v / A_net  # Velocity of vapour flowing up the column
        Flv = L_w / V_w * (rho_v / rho_l) ** 0.5  # Eq 11.82
        u_f = table_K1(Flv, l_spacing) * ((rho_l - rho_v) / rho_v) ** 0.5 * (sigma / 0.02) ** 0.2  # * min(1, 5*A_hs/A_active)

        flooding = u_f - u

        # entrainment
        percent_flooding = u / u_f * 100
        psi = table_psi(Flv, percent_flooding)

        entrainment = 0.1 - psi

        # weeping
        h_ow = 0.75 * (L_w / rho_l / l_weir) ** (2 / 3)
        Z = h_ow + h_weir
        K2 = table_K2(1000 * Z)
        u_h_min = (K2 - 0.9 * (25.4 - d_hole * 1000)) / rho_v ** 0.5
        u_h = V_w / rho_v / A_hs

        weeping = u_h - u_h_min

        # pressure calculations
        C0 = table_C0(A_hs / A_p, t_plate / d_hole)
        h_d = 0.051 * rho_v / rho_l * (u_h / C0) ** 2  # m
        h_r = 12.5 / rho_l  # m
        h_t = h_d + h_ow + h_weir + h_r

        # down comer backup
        A_m = min(A_down, l_weir * h_ap)  # ideally would not have this so as to keep function smooth
        h_dc = 0.166 * (L_w / rho_l / A_m) ** 2  # m
        h_b = h_weir + h_ow + h_t + h_dc

        backup = 0.5 * (h_weir + l_spacing) - h_b

        # residence time
        t_r = A_down * h_b * rho_l / L_w

        residence = t_r - 3
        return flooding, entrainment, weeping, backup, residence
    
    return _constraints(V, L) + _constraints(V*turn_down, L*turn_down)

def check_design(parameters, design):
    names = CONSTRAINT_NAMES + [f'{n} at turndown' for n in CONSTRAINT_NAMES]
    values = columnconstraints(parameters, design)
    for v, n in zip(values, names):
        status = '✅' if v > 0 else '🚫'
        print(f'{status} {n}: {v}')

if __name__ == "__main__":
    V_w = 55.5 * 55.6 / 3600
    L_w = 811.6 * 18 / 3600
    rho_v = 0.72
    rho_l = 954
    mu = 1e-4  # could not find value in C&R, not used in constraints though
    sigma = 57e-3
    turn_down = 0.7
    t_plate = 5e-3
    l_calm = 50e-3

    d_col = 0.79
    d_hole = 5e-3
    theta = 99 / 180 * numpy.pi
    l_spacing = 0.5
    l_pitch = 12.5e-3
    h_weir = 50e-3
    h_ap = 40e-3

    parameters = [V_w, L_w, rho_v, rho_l, mu, sigma, turn_down, t_plate, l_calm]
    design = [d_col, d_hole, theta, l_spacing, l_pitch, h_weir, h_ap]

    check_design(parameters, design)

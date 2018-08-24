from pyomo.environ import sin, cos, acos, sqrt, log, log10, exp
from numpy import pi
from costs import materials
from classes import Parameters, Design


def table_K1(Flv, l_spacing):
    """Figure 11-27
    
    Regression by Dazzles
    """
    a, b, c, d, e, f, g = -0.4141022057164735, -0.15097503976930754, -0.03955953480149698, -0.6833440064263923, 1.2631972777123914,       -0.10683783034172412, -1.8165949194259345
    return 10**(a*log10(Flv) + b*(log10(Flv))**2 + c*exp(Flv) + d*l_spacing**2 + e*l_spacing + f*log10(Flv)*l_spacing + g)


def table_psi(Flv, flood):
    """Figure 11-29

    Regression by Dazzles
    """
    a, b, c, d, e, f = -1.7890463988780416, 0.3276578078888633, -0.00016483063266354298, 0.029666898603030227, -0.011235456340583339, -2.8085647498602717
    return 10**(a*exp(log10(Flv)) + b*log10(Flv) + c*flood**2 + d*flood + e*log10(Flv)*flood + f)


def table_K2(h):
    """Figure 11-30
    
    Regression by Carl Sandrock"""
    a, b, c, d = 26.52258403,  0.76197575, 13.23116619, 33.1867269
    return a + b*sqrt(sqrt((h - c)**2)) - (h - c)/d


def table_C0(A_ratio, P_ratio):
    """Figure 11-34

    Regression by Dazzles
    """
    a, b, c, d = 0.007910348009503125, 0.1618154799051229, -0.0395502040700688, 0.6334016756405921
    return a*A_ratio + b*P_ratio**2 + c*P_ratio + d


CONSTRAINT_NAMES = ["Flooding",
                    "Entrainment",
                    "Weeping",
                    "Downcomer backup",
                    "Residence time"]


def columnconstraints(parameters, design):
    L, V, rho_l, rho_v, sigma, turn_down, t_shell, l_calm, Nplates = parameters
    d_col, theta, d_hole, l_pitch, l_spacing, h_weir, h_ap, t_plate = design

    # Geometry calculations
    l_weir = d_col*sin(theta/2)  # calculated using trig

    h = d_col/2*cos(theta/2)
    h1 = h - l_calm
    d1 = d_col - 2*l_calm
    
    theta1 = 2*cos(h1/(d1/2))
    A1 = (theta1 - sin(theta1))/2*(d1/2)**2

    # Areas
    A_col = pi/4*d_col**2
    A_h = pi/4*d_hole**2
    A_down = (theta - sin(theta))/2*(d_col/2)**2
    # A_active = A_col - 2 * A_down
    A_p = pi/4*d1**2 - 2*A1
    A_hs = (A_h / 2)/(sqrt(3)/4*l_pitch**2)*A_p
    A_net = A_col - A_down

    def _constraints(V_w, L_w):
        # Flooding
        u = V_w/rho_v/A_net  # Velocity of vapour flowing up the column
        Flv = L_w/V_w*(rho_v/rho_l) ** 0.5  # Eq 11.82
        u_f = table_K1(Flv, l_spacing) * ((rho_l - rho_v)/rho_v)**0.5*(sigma/0.02)**0.2  # * min(1, 5*A_hs/A_active)

        flooding = u_f - u

        # entrainment
        percent_flooding = u/u_f*100
        psi = table_psi(Flv, percent_flooding)

        entrainment = 0.1 - psi

        # weeping
        h_ow = 0.75*(L_w/rho_l/l_weir)**(2/3)
        Z = h_ow + h_weir
        K2 = table_K2(1000*Z)
        u_h_min = (K2 - 0.9*(25.4 - d_hole * 1000))/rho_v**0.5
        u_h = V_w/rho_v/A_hs

        weeping = u_h - u_h_min

        # pressure calculations
        C0 = table_C0(A_hs/A_p, t_plate/d_hole)
        h_d = 0.051*rho_v/rho_l*(u_h/C0)**2  # m
        h_r = 12.5/rho_l  # m
        h_t = h_d + h_ow + h_weir + h_r

        # down comer backup
        A_m = l_weir * h_ap
        h_dc = 0.166*(L_w/rho_l/A_m)**2  # m
        h_b = h_weir + h_ow + h_t + h_dc

        backup = 0.5*(h_weir + l_spacing) - h_b

        # residence time
        t_r = A_down*h_b*rho_l/L_w

        residence = t_r - 3
        return flooding, entrainment, weeping, backup, residence
    
    return _constraints(V, L) + _constraints(V*turn_down, L*turn_down)


def check_design(parameters, design):

    print("Parameters:")
    print(parameters)

    print("Design:")
    print(design)

    names = CONSTRAINT_NAMES + [f'{n} at turndown' for n in CONSTRAINT_NAMES]
    values = columnconstraints(parameters, design)

    print("Constraints:")
    for v, n in zip(values, names):
        status = '✅' if v > 0 else '🚫'
        print(f'{status} {n}: {v}')

    print("Cost:")

    print(materials(parameters, design))


if __name__ == "__main__":
    # We check the design in Example 11.2 and Example 11.11 in C+R

    # Note in the design they design the bottom plate only, but we assume
    # the same plate design throughout
    parameters = Parameters(
        L_w=811.6*18/3600,
        V_w=55.5*55.6/3600,
        rho_l=954,
        rho_v=0.72,
        sigma=57e-3,
        turn_down=0.7,
        t_shell=6e-3,
        l_calm=50e-3,
        Nplates=25,
    )

    design = Design(
        d_col=0.79,
        theta=99/180*pi,
        d_hole=5e-3,
        l_pitch=12.5e-3,
        l_spacing=0.5,
        h_weir=50e-3,
        h_ap=40e-3,
        t_plate=5e-3,
    )

    check_design(parameters, design)

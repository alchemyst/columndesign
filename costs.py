import numpy

def materials(parameters, design):

    diameter = design.d_col
    shell_thickness = parameters.t_shell
    plate_thickness = design.t_plate
    angle = design.theta
    plates = parameters.Nplates
    spacing = design.l_spacing

    r = diameter/2
    shell_cost = (numpy.pi*diameter)*(plates*spacing)*shell_thickness  # circumference * height * thickness
    plate_area = numpy.pi*r**2
    downcomer_area = (angle - numpy.sin(angle))*r**2/2
    plate_cost = plates*plate_thickness*(plate_area - downcomer_area)

    return shell_cost + plate_cost

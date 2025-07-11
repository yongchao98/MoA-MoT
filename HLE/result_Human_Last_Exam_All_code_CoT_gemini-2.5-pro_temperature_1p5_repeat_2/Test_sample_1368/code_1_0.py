def solve_construction():
    """
    This function determines and prints the shortest sequence of commands
    to construct a square inscribed in a circle.
    
    L: Draw a line between two points.
    C: Draw a circle with a center and a point on its radius.
    
    The steps are:
    1. L: Draw a line through the center and the given point on the circumference to create a diameter. This defines the first two vertices.
    2. C: Draw a circle centered on the first vertex with a radius equal to the diameter.
    3. C: Draw a second circle centered on the second vertex with the same diameter radius.
    4. L: Draw a line between the intersection points of the two new circles. This line is the perpendicular bisector of the diameter, creating a second diameter.
    
    The four points where the two diameters intersect the original circle are the vertices of the inscribed square.
    """
    
    # The shortest sequence of commands
    command_sequence = "LCCL"
    
    print(command_sequence)

solve_construction()
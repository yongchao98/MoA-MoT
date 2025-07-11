import math

def solve_rod_sliding_angle():
    """
    This function prints the derived formula for the angle at which the rod begins to slide.
    The formula is symbolic, expressed in terms of the given parameters.
    """
    
    # The parameters are:
    # L: total length of the rod
    # M: total mass of the rod
    # l: distance from the pivot (table edge) to the center of mass of the rod
    # mu: coefficient of static friction
    # theta: the angle of the rod with the horizontal
    
    print("The derivation shows that sliding occurs when the friction force 'f' equals 'mu * N'.")
    print("The forces depend on the dynamics of the rotation.")
    print("\nThe relationship is given by tan(theta) = mu * FACTOR.")
    print("The derived expression for tan(theta) is:")
    
    # We construct the string for the expression part by part to meet the output requirements.
    # The numerator of the factor multiplying mu
    numerator = "(L**2 / 12) + (2 * l**2)"
    
    # The denominator of the factor multiplying mu
    denominator = "(L**2 / 12) + (3 * l**2)"
    
    # Print the equation for tan(theta)
    print(f"\ntan(theta) = mu * ( {numerator} ) / ( {denominator} )\n")
    
    # Print the final expression for the angle theta
    print("Therefore, the angle theta at which the rod begins to slide is:")
    print(f"\ntheta = arctan( mu * ( {numerator} ) / ( {denominator} ) )\n")

solve_rod_sliding_angle()

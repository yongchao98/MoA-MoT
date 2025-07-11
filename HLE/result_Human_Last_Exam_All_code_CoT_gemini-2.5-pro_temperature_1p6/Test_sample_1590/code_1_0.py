import math

def print_angle_expression():
    """
    Prints the symbolic expression for the angle theta at which the rod begins to slide.
    
    The variables in the expression are:
    L: the length of the rod
    l: the distance from the pivot (table edge) to the center of mass of the rod
    mu: the coefficient of static friction between the rod and the table edge
    """
    
    # The final derived expression for tan(theta) is:
    # tan(theta) = (mu * L**2 + 12 * l**2 * (1 + mu)) / (L**2 + 12 * l**2 * (1 + 2 * mu))
    # The expression for theta is the arctan of this expression.
    
    final_expression = "theta = arctan( (mu * L**2 + 12 * l**2 * (1 + mu)) / (L**2 + 12 * l**2 * (1 + 2 * mu)) )"
    
    print("The expression for the angle theta at which the rod begins to slide is:")
    print(final_expression)

if __name__ == "__main__":
    print_angle_expression()
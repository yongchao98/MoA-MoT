import math

def display_period_formula():
    """
    Prints the derived formula for the period of motion of the system.
    """
    # These are the numerical coefficients derived from the Lagrangian analysis.
    coefficient_outside_sqrt = 2
    numerator_R_coefficient = 101
    denominator_g_coefficient = 40

    print("The period of motion, T, is given by the following equation:")
    print("")
    # Using f-string to format the output clearly
    print(f"T = {coefficient_outside_sqrt} * π * sqrt( ({numerator_R_coefficient} * R) / ({denominator_g_coefficient} * g) )")
    print("")
    print("Where:")
    print("  π  is the mathematical constant pi (approx. 3.14159)")
    print("  R  is the radius of the identical disks")
    print("  g  is the acceleration due to gravity")

# Execute the function to display the result
display_period_formula()
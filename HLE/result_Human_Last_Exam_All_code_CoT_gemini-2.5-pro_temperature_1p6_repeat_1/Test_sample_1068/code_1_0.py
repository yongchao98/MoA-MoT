import math

def solve_l_14():
    """
    This function calculates the value of l(14) based on the derived formula
    l(14) = 2 * ln( Gamma(0.5) * Gamma(15) / (Gamma(7.5) * Gamma(8)) )
    and outputs the components of the final simplified equation.
    """

    # Parameters for the Gamma function formula derived from the integral
    gamma_arg_1 = 0.5
    gamma_arg_2 = 15.0
    gamma_arg_3 = 7.5
    gamma_arg_4 = 8.0

    # Calculate the values of the Gamma functions
    val_gamma_1 = math.gamma(gamma_arg_1)
    val_gamma_2 = math.gamma(gamma_arg_2)
    val_gamma_3 = math.gamma(gamma_arg_3)
    val_gamma_4 = math.gamma(gamma_arg_4)
    
    # The final equation is derived as l(14) = 2 * ln(P)
    # where P = Gamma(0.5) * Gamma(15) / (Gamma(7.5) * Gamma(8))
    # We first print the numbers for this equation.
    print(f"The equation for l(14) in terms of the Gamma function is:")
    print(f"l(14) = 2 * ln( Gamma({gamma_arg_1}) * Gamma({gamma_arg_2}) / (Gamma({gamma_arg_3}) * Gamma({gamma_arg_4})) )")
    print("\nCalculating the values:")
    print(f"Gamma({gamma_arg_1}) = {val_gamma_1}")
    print(f"Gamma({gamma_arg_2}) = {val_gamma_2}")
    print(f"Gamma({gamma_arg_3}) = {val_gamma_3}")
    print(f"Gamma({gamma_arg_4}) = {val_gamma_4}")

    # Calculate the value of the expression P inside the logarithm
    P = (val_gamma_1 * val_gamma_2) / (val_gamma_3 * val_gamma_4)
    
    print(f"\nThe value of the expression inside the logarithm is: {P}")

    # The expression simplifies to 28 * ln(2)
    # The numbers in this final equation are 28 and 2.
    final_coeff_1 = 28
    final_coeff_2 = 2

    print(f"\nThe expression simplifies to the final equation: {final_coeff_1} * ln({final_coeff_2})")
    
    # Calculate the final numerical value
    final_value = final_coeff_1 * math.log(final_coeff_2)
    print(f"\nThe numerical value of l(14) is: {final_value}")

solve_l_14()
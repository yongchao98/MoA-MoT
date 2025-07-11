import math

def final_answer():
    """
    This function prints the normalised density of the invariant measure for the given map.
    """
    
    # The functional form of the density is rho(x) = C / (x + x^(1/2))
    # The normalization constant C is 1 / (2 * ln(2))
    
    coeff_of_ln = 2
    arg_of_ln = 2
    
    # The powers of x in the denominator
    power_1 = 1.0
    power_2 = 0.5
    
    constant_value = coeff_of_ln * math.log(arg_of_ln)
    
    print("The normalised density of the invariant measure is given by the formula:")
    print(f"rho(x) = 1 / (({coeff_of_ln} * ln({arg_of_ln})) * (x^{power_1} + x^{power_2}))")
    print("\nWhere the numerical constants in the formula are:")
    print(f"Coefficient of the logarithm: {coeff_of_ln}")
    print(f"Argument of the logarithm: {arg_of_ln}")
    print(f"Powers of x in the denominator: {power_1} and {power_2}")
    
    print(f"\nThe value of the normalization constant {coeff_of_ln}*ln({arg_of_ln}) is approximately {constant_value:.6f}")

final_answer()
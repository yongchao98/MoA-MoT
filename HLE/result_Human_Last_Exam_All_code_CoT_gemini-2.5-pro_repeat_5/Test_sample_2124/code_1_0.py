import sympy

def calculate_R_ratio():
    """
    Calculates the ratio R based on the derived coefficients of the counter-terms.
    
    The counter-terms are proportional to a common factor C = g^2 / (32 * pi^2 * epsilon).
    dZ_x = 1 * C
    dZ_mx = -2 * C
    dZ_g = -3 * C
    
    The ratio R is independent of this common factor.
    """
    
    # Coefficients of the common factor C for each counter-term
    dZ_x_coeff = 1
    dZ_mx_coeff = -2
    dZ_g_coeff = -3
    
    # Calculate the numerator and denominator
    numerator = dZ_x_coeff
    denominator = dZ_g_coeff + dZ_mx_coeff
    
    # Calculate the ratio R
    R = sympy.Rational(numerator, denominator)
    
    # Print the equation with the numbers
    print(f"The ratio R is calculated as:")
    print(f"R = dZ_x / (dZ_g + dZ_mx)")
    # Representing the coefficients in the equation
    print(f"R = ({dZ_x_coeff}) / (({dZ_g_coeff}) + ({dZ_mx_coeff}))")
    print(f"R = {numerator} / {denominator}")
    print(f"R = {R}")

calculate_R_ratio()
import math
from fractions import Fraction

def get_p_formula():
    """
    Calculates the coefficients for the P(n) formula based on the Euler-Maclaurin expansion
    and prints the resulting formula for P(n).
    """

    # Bernoulli numbers needed for the expansion
    B4 = Fraction(-1, 30)
    B6 = Fraction(1, 42)

    # Derivatives of f(x) = x*ln(x) evaluated at n as symbolic expressions
    # f'''(n) = -1/n^2
    # f^(5)(n) = -6/n^4
    # We only need the constant factors of these derivatives.
    f3_coeff = -1
    f5_coeff = -6

    # Calculate the coefficients for the terms in the log expansion
    # Term for n^-2: (B4/4!) * f'''(n)
    coeff1 = (B4 / math.factorial(4)) * f3_coeff
    # Term for n^-4: (B6/6!) * f^(5)(n)
    coeff2 = (B6 / math.factorial(6)) * f5_coeff
    
    # Denominators for the final formula
    den1 = coeff1.denominator
    den2 = coeff2.denominator
    
    # We expect coeff2 to be negative, so we print its absolute value after a minus sign.
    print(f"P(n) = exp((1/{den1})/n^2 - (1/{abs(den2)})/n^4)")
    
    # We print the numbers in the equation for clarity as requested.
    print("\nCalculated numbers for the equation:")
    print(f"Coefficient 1 numerator: {coeff1.numerator}")
    print(f"Coefficient 1 denominator: {den1}")
    print(f"Coefficient 2 numerator: {coeff2.numerator}")
    print(f"Coefficient 2 denominator: {den2}")

if __name__ == '__main__':
    get_p_formula()

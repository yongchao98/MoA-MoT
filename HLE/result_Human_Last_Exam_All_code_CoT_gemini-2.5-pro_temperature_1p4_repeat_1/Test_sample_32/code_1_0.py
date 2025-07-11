import math
from fractions import Fraction

def get_coeffs(n):
    """
    Computes the coefficients C_j = (6j)! / ((3j)! * (2j)!) for j from 0 to n.
    """
    coeffs = []
    for j in range(n + 1):
        # Use logs for stability with large numbers before converting to Fraction
        if j == 0:
            coeffs.append(Fraction(1))
            continue
        # Numerator is math.factorial(6*j)
        # Denominator is math.factorial(3*j) * math.factorial(2*j)
        # We can compute this with integer arithmetic using Fraction for precision
        c = Fraction(math.factorial(6*j), math.factorial(3*j) * math.factorial(2*j))
        coeffs.append(c)
    return coeffs

def compute_integral():
    """
    Computes the integral of lambda_3*lambda_2*lambda_1 on M_3.
    """
    g = 3
    degree = 3 * g - 3 # We need the coefficient of z^6
    
    # We need coefficients C_j for j from 0 to 6
    C = get_coeffs(degree)
    
    total_sum = Fraction(0)
    
    # Sum over all j, k, l such that j + k + l = 6
    for j in range(degree + 1):
        for k in range(degree - j + 1):
            l = degree - j - k
            
            term_value = C[j] * C[k] * C[l]
            
            # The powers are from (1*z), (2*z), (3*z)
            # The powers correspond to lambda_1, lambda_2, lambda_3
            term_value *= Fraction(1**j * 2**k * 3**l)
            
            total_sum += term_value

    # The formula has a prefactor of (1/36)^6
    final_integral = total_sum * Fraction(1, 36**degree)
    
    return final_integral

result = compute_integral()

# The final result is presented in the format "a/b"
# The numerator and denominator of the fraction are printed in the final equation.
print(f"The integral of lambda_3 * lambda_2 * lambda_1 on the moduli space of genus 3 curves is:")
print(f"Integral = {result.numerator}/{result.denominator}")
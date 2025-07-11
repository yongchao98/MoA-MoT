import numpy as np
from scipy.integrate import quad
from fractions import Fraction

def solve_integral():
    """
    This function calculates the specified integral numerically and returns the
    result as a fraction.
    """

    # We define the original, unsimplified integrand to be robust against
    # potential derivation errors.
    def integrand(x):
        # A small epsilon to handle potential floating point issues at boundaries
        # although for this problem it's not strictly necessary.
        eps = 1e-9
        
        # Ensure x is not exactly 0 or pi to avoid warnings, though abs handles it.
        x_safe = np.clip(x, eps, np.pi - eps)

        sin_x = np.sin(x_safe)
        cos_2x = np.cos(2 * x_safe)
        sin_2x = np.sin(2 * x_safe)
        cos_3x = np.cos(3 * x_safe)
        
        term_A = np.abs(2 * sin_x)
        term_B = np.abs(2 * cos_2x - 1)
        
        term_C = np.abs(sin_2x)
        term_D = np.abs(cos_3x)
        
        max_val = np.maximum(term_A, term_B)
        min_val = np.minimum(term_C, term_D)
        
        return (max_val * min_val)**50

    # Perform the numerical integration from 0 to pi.
    # We increase the limit and decrease the tolerance for high precision.
    integral_value, error = quad(integrand, 0, np.pi, limit=200, epsabs=1e-15, epsrel=1e-15)

    # The result is expected to be a rational multiple of pi.
    # We find the coefficient by dividing by pi.
    coefficient = integral_value / np.pi

    # Use the Fraction class to find the exact rational number.
    # We set a large denominator limit to capture the precise fraction.
    frac = Fraction(coefficient).limit_denominator(10**18)
    
    numerator = frac.numerator
    denominator = frac.denominator
    
    # The final equation is I = (numerator / denominator) * pi
    # The user asks for the answer "strictly as a fraction", which we interpret
    # as the coefficient of pi.
    print("The integral evaluates to an expression of the form: (p/q) * pi")
    print(f"The fraction p/q has been calculated as:")
    print(f"p = {numerator}")
    print(f"q = {denominator}")
    print("\nThe final equation is:")
    print(f"I = ( {numerator} / {denominator} ) * pi")
    
    # The required answer format is the fraction string itself.
    final_answer_fraction = f"{numerator}/{denominator}"
    
    return final_answer_fraction

if __name__ == '__main__':
    answer = solve_integral()
    print("\nAnswer in the requested format:")
    print(f"<<<{answer}>>>")

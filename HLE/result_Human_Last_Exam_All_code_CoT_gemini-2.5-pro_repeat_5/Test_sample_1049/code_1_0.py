import sys
from fractions import Fraction

def get_closed_form_string():
    """
    This function returns the string representation of the closed form for the sum.
    The coefficients are derived mathematically as explained above.
    """
    coeffs = [
        Fraction(63, 8),
        Fraction(245, 8),
        Fraction(355, 8),
        Fraction(235, 8),
        Fraction(70, 8),
        Fraction(8, 8)
    ]

    # Build the polynomial string
    poly_terms = []
    for i, c in enumerate(coeffs):
        power = 5 - i
        # Format the coefficient as a fraction if it has a denominator other than 1
        if c.denominator == 1:
            coeff_str = f"{c.numerator}"
        else:
            coeff_str = f"({c})"

        # Format the power of n
        if power > 1:
            term = f"{coeff_str} * n^{power}"
        elif power == 1:
            term = f"{coeff_str} * n"
        else:
            term = f"{coeff_str}"
        
        poly_terms.append(term)
    
    poly_str = " + ".join(poly_terms)
    
    return f"4^n * ({poly_str})"

if __name__ == '__main__':
    # The final closed form is an equation in terms of n.
    # The code prints this equation.
    final_equation = get_closed_form_string()
    print("The closed form for the sum is:")
    print(final_equation)

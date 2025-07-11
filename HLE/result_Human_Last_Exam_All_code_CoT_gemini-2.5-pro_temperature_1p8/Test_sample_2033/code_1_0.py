import math

def calculate_l_value(a, b, c, d):
    """
    This function calculates the value of l(a,b,c,d) based on the derived
    analytical formula.
    """
    n = 20

    def S(x, num_terms):
        """Calculates the sum of a geometric series: x^1 + ... + x^num_terms."""
        if abs(x - 1.0) < 1e-12:
            return num_terms
        if abs(1.0 - x) < 1e-12:
            return float('nan') # Should not happen with problem constraints
        return x * (1 - x**num_terms) / (1 - x)

    # The formula for l is derived from 0.5 * Tr(K_a^2 * D_b^-1 * (D_d - D_c))
    # After expanding and summing the geometric series, we get:
    # l = 1/(2*(1-a^2)) * [ term1 - term2 - term3 ]

    # Term 1: from (1+a^2) part of diag(K_a^2)
    s_db = S(d / b, n)
    s_cb = S(c / b, n)
    term1 = (1 + a**2) * (s_db - s_cb)

    # Term 2: from -a^(2i) part of diag(K_a^2)
    s_a2db = S(a**2 * d / b, n)
    s_a2cb = S(a**2 * c / b, n)
    term2 = s_a2db - s_a2cb

    # Term 3: from -a^(2(n+1-i)) part of diag(K_a^2)
    s_a2bd = S(a**2 * b / d, n)
    s_a2bc = S(a**2 * b / c, n)
    term3 = ((d / b)**(n + 1)) * s_a2bd - ((c / b)**(n + 1)) * s_a2bc
    
    # Final combination
    numerator = term1 - term2 - term3
    denominator = 2 * (1 - a**2)
    
    l_val = numerator / denominator

    # Printing the components of the final equation as requested
    print("This calculation is based on the simplified analytical formula:")
    print("l = (Term1 - Term2 - Term3) / (2 * (1 - a^2))\n")
    print(f"For parameters a={a}, b={b}, c={c}, d={d}, n={n}:\n")
    print("--- Components of the Final Equation ---")
    print(f"Term1 = (1+a^2) * (S(d/b, n) - S(c/b, n))                     = {term1}")
    print(f"Term2 = S(a^2*d/b, n) - S(a^2*c/b, n)                          = {term2}")
    print(f"Term3 = (d/b)^(n+1)*S(a^2*b/d, n) - (c/b)^(n+1)*S(a^2*b/c, n) = {term3}")
    print(f"Numerator (Term1 - Term2 - Term3)                            = {numerator}")
    print(f"Denominator (2 * (1 - a^2))                                    = {denominator}")
    print("----------------------------------------\n")

    print("Final Equation:")
    print(f"l = ({term1:.4f} - {term2:.4f} - {term3:.4f}) / {denominator:.4f}")
    print(f"l = {numerator:.4f} / {denominator:.4f}\n")
    
    print("Result:")
    print(l_val)
    

# Example usage with some values satisfying the constraints.
# a in (0, 1), b,c,d in [1, inf]
a_in = 0.5
b_in = 1.5
c_in = 2.0
d_in = 2.5

calculate_l_value(a_in, b_in, c_in, d_in)

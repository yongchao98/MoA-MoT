from fractions import Fraction

def solve():
    """
    This function calculates the coefficients for the correction factor P(n)
    and prints the final formula.
    """
    # The asymptotic expansion for ln(Q(n)/T(n)) is given by:
    # ln(Q(n)/T(n)) ~ c_2*n^-2 + c_4*n^-4 + O(n^-6)
    # where c_2 = -1/720 and c_4 = 1/5040.
    
    # We want to find P(n) such that ln(Q(n)/(T(n)P(n))) = O(n^-6).
    # This requires ln(P(n)) to match the expansion of ln(Q(n)/T(n)) up to the relevant terms.
    # ln(P(n)) = c_2*n^-2 + c_4*n^-4
    ln_p_c2 = Fraction(-1, 720)
    ln_p_c4 = Fraction(1, 5040)
    
    # We find P(n) by expanding P(n) = exp(ln(P(n))).
    # P(n) = 1 + (c_2*n^-2 + c_4*n^-4) + 1/2 * (c_2*n^-2 + c_4*n^-4)^2 + ...
    # We need to collect terms for the series of P(n) = 1 + a_2*n^-2 + a_4*n^-4 + ...
    
    # Coefficient for n^-2
    a2 = ln_p_c2
    
    # Coefficient for n^-4 is ln_p_c4 + (ln_p_c2)^2 / 2
    a4 = ln_p_c4 + (ln_p_c2**2) / 2
    
    # The terms with odd powers (n^-1, n^-3, n^-5) are zero.
    # So P(n) is truncated after the n^-4 term to achieve O(n^-6) error.
    
    # The coefficients are:
    a2_num = a2.numerator
    a2_den = a2.denominator
    
    a4_num = a4.numerator
    a4_den = a4.denominator
    
    # The final equation string needs to show each number.
    # We handle the sign of the coefficients explicitly.
    sign_a2 = "-" if a2_num < 0 else "+"
    abs_a2_num = abs(a2_num)
    
    sign_a4 = "-" if a4_num < 0 else "+"
    abs_a4_num = abs(a4_num)

    # Note that the problem implies we should format the formula as requested.
    # The final P(n) is the expansion up to n^-5. Since there are only even powers, this includes the n^-4 term.
    
    print("The formula for P(n) is:")
    print(f"P(n) = 1 {sign_a2} ({abs_a2_num}/{a2_den})/n^2 {sign_a4} ({abs_a4_num}/{a4_den})/n^4")

solve()

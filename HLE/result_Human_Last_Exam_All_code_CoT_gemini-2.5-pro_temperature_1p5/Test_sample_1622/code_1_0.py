def solve_p_n_formula():
    """
    This function derives and prints the formula for P(n) based on asymptotic analysis.

    The problem is to find P(n) such that the approximation:
    T_new(n) = A * n**( (ln(n))/2 ) * (1 + ln(n)/(2*n) + P(n))
    has a relative error of O(((ln n)/n)**4) for Q(n).
    Let L = ln(n).

    This requires matching the asymptotic expansion of ln(Q(n)) with that of ln(T_new(n)).

    The expansion of ln(Q(n)) - (ln(A) + L**2/2) is:
    L/(2*n) + (1-L)/(12*n**2) + 0/n**3 + O(L/n**4)

    The expansion of ln(1 + L/(2*n) + P(n)) is matched against this.
    Let P(n) = C2/n**2 + C3/n**3 + ...

    Matching the 1/n**2 term gives C2:
    C2 - L**2/8 = (1-L)/12
    => C2 = (3*L**2 - 2*L + 2) / 24

    Matching the 1/n**3 term gives C3:
    C3 - L*C2/2 + L**3/24 = 0
    => C3 = L*C2/2 - L**3/24 = (L**3 - 2*L**2 + 2*L) / 48

    These terms for P(n) are sufficient to make the error O((L/n)**4).
    The code below constructs the formula string from these derived coefficients.
    """

    # Coefficients for the n^-2 term numerator: 3L^2 - 2L + 2
    p2_coeffs = {'L**2': 3, 'L': -2, 'const': 2}
    p2_denom = 24

    # Coefficients for the n^-3 term numerator: L^3 - 2L^2 + 2L
    p3_coeffs = {'L**3': 1, 'L**2': -2, 'L': 2}
    p3_denom = 48

    # Build the string for the first term of P(n)
    p2_num_str = f"({p2_coeffs['L**2']}*L**2 {p2_coeffs['L']:+}*L + {p2_coeffs['const']})"
    p2_term_str = f"{p2_num_str}/({p2_denom}*n**2)"

    # Build the string for the second term of P(n)
    p3_num_str = f"({p3_coeffs['L**3']}*L**3 {p3_coeffs['L**2']:+}*L**2 + {p3_coeffs['L']}*L)"
    p3_term_str = f"{p3_num_str}/({p3_denom}*n**3)"

    # Combine the terms for the full P(n) formula
    formula = f"P(n) = {p2_term_str} + {p3_term_str}"

    # Clean up the string for better readability
    formula = formula.replace('+ -', '- ')
    formula = formula.replace('(1*L', '(L')

    print(formula)

solve_p_n_formula()
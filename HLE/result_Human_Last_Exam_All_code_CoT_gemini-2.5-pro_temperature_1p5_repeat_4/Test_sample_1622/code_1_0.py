def solve():
    """
    This function prints the formula for P(n).
    L is used as a placeholder for ln(n).
    """
    # Coefficients for the n^-2 term: (a*L^2 + b*L + c)
    a2 = 3
    b2 = -2
    c2 = 2
    den2 = 24

    # Coefficients for the n^-3 term: (a*L^3 + b*L^2 + c*L)
    a3 = 1
    b3 = -2
    c3 = 2
    den3 = 48

    # Using f-string to format the mathematical formula for P(n)
    # The numbers are explicitly written as requested by the prompt.
    formula = f"P(n) = ({a2}*L^2 - {abs(b2)}*L + {c2}) / ({den2}*n^2) + ({a3}*L^3 - {abs(b3)}*L^2 + {c3}*L) / ({den3}*n^3)"

    print(formula)

solve()
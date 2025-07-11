def solve():
    """
    This function prints the formula for P(n) which refines the approximation for Q(n).
    L is used as a shorthand for ln(n).
    """
    
    # Coefficients for the 1/n^2 term
    c22_num = 3
    c21_num = 2
    c20_num = -2
    den2 = 24
    
    # Coefficients for the 1/n^3 term
    c33_num = 1
    c32_num = 2
    c31_num = -2
    den3 = 48

    # Constructing the formula string
    # The instruction is to output each number in the final equation.
    # So we build the string with the explicit numbers.
    
    term2 = f"({c22_num}L^2 + {c21_num}L {c20_num})/({den2}n^2)"
    term3 = f"({c33_num}L^3 + {c32_num}L^2 {c31_num}L)/({den3}n^3)"
    
    formula = f"P(n) = {term2} + {term3}"
    
    print(formula)

solve()
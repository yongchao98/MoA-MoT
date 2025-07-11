import sympy

def derive_P_n():
    """
    This function derives and prints the formula for P(n).
    L is a symbol for ln(n).
    """
    n = sympy.Symbol('n')
    L = sympy.Symbol('L')

    # According to the derivation explained in the plan, P(n) is found to be a sum
    # of terms in powers of 1/n. To achieve the required accuracy, we need the terms
    # up to n^(-3).

    # The coefficient of the n^-2 term is (3*L^2 - 2*L + 2) / 24.
    c2_num = [3, -2, 2]
    c2_den = 24

    # The coefficient of the n^-3 term is (L^3 - 2*L^2 + 2*L) / 48.
    c3_num = [1, -2, 2, 0] # Coefficients for L^3, L^2, L^1, L^0
    c3_den = 48
    
    print("The formula for P(n) is derived as P(n) = P_2(n) + P_3(n) where L = ln(n).")
    print("\nThe term P_2(n) is:")
    print(f"P_2(n) = (({c2_num[0]})*L**2 + ({c2_num[1]})*L + ({c2_num[2]})) / ({c2_den}*n**2)")

    print("\nThe term P_3(n) is:")
    print(f"P_3(n) = (({c3_num[0]})*L**3 + ({c3_num[1]})*L**2 + ({c3_num[2]})*L) / ({c3_den}*n**3)")
    
    print("\nCombining these gives the full formula for P(n):")
    # Using sympy to pretty print the final formula.
    p2_expr = (c2_num[0]*L**2 + c2_num[1]*L + c2_num[2]) / (c2_den * n**2)
    p3_expr = (c3_num[0]*L**3 + c3_num[1]*L**2 + c3_num[2]*L) / (c3_den * n**3)
    P_n_expr = p2_expr + p3_expr
    
    sympy.pprint(sympy.Eq(sympy.Symbol("P(n)"), P_n_expr), use_unicode=False)

if __name__ == '__main__':
    derive_P_n()
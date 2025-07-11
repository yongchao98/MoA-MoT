def solve_svm_bound_coeffs():
    """
    This function determines and prints the coefficients c1 and c2 for the extended Jaakola-Haussler bound.
    """
    # Based on the derivation, the coefficients are:
    c1 = 1
    c2 = 0

    print("The determined coefficients are:")
    print(f"c1 = {c1}")
    print(f"c2 = {c2}")
    
    print("\nThe final inequality is:")
    print(f"-(K * alpha_D_minus_i)_i <= (1 + {c1}*beta)*alpha_D_i - (1 + {c2}*beta)*(K*alpha_D)_i + o(beta)")
    print("Which simplifies to:")
    print(f"-(K * alpha_D_minus_i)_i <= (1 + beta)*alpha_D_i - (K*alpha_D)_i + o(beta)")

solve_svm_bound_coeffs()
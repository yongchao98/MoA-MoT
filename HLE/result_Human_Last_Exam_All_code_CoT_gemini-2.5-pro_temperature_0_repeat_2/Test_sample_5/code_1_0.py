def print_proportionality_factor_formula():
    """
    This function prints the coefficients of the expanded formula for the proportionality factor C(d, k).
    The formula is C(d, k) = d - (d - 2k)^2 = -d^2 + (4k+1)d - 4k^2.
    """
    
    # Coefficients of the expanded formula: c_d2 * d^2 + c_d * d + c_dk * dk + c_k2 * k^2
    c_d2 = -1
    c_d = 1
    c_dk = 4
    c_k2 = -4
    
    print("The proportionality factor is C(d, k) = d - (d - 2k)^2.")
    print("Expanding this gives: C(d, k) = c_d*d + c_d2*d^2 + c_dk*d*k + c_k2*k^2")
    print("The coefficients are:")
    print(f"Coefficient of d: {c_d}")
    print(f"Coefficient of d^2: {c_d2}")
    print(f"Coefficient of d*k: {c_dk}")
    print(f"Coefficient of k^2: {c_k2}")
    
    # The prompt asks to output each number in the final equation.
    # Let's print the expanded equation with its numerical coefficients.
    print("\nThe final equation with numerical coefficients is:")
    print(f"C(d, k) = ({c_d})*d + ({c_d2})*d^2 + ({c_dk})*d*k + ({c_k2})*k^2")

print_proportionality_factor_formula()
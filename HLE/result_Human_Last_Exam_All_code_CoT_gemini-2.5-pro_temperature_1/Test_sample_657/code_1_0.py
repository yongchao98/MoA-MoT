import sympy

def generate_controller_expression():
    """
    This function generates and prints the mathematical expression for the set of
    all proper stabilizing controllers H_2(s) for the plant H_1(s) = s/(s^2-1).
    The controller is parametrized by K(s), which can be any stable and proper
    transfer function.
    """
    # Define the symbolic variables
    s = sympy.Symbol('s')
    K = sympy.Function('K')(s)

    # Numerator of the controller H_2(s)
    # This comes from (s^2 - 1) * K(s)
    num_poly_part = sympy.poly(s**2 - 1, s)
    num_coeffs = num_poly_part.all_coeffs() # [1, 0, -1]

    # Denominator of the controller H_2(s)
    # This comes from (s^2 + 2s + 1) - s*K(s)
    den_poly_part = sympy.poly(s**2 + 2*s + 1, s)
    den_coeffs = den_poly_part.all_coeffs() # [1, 2, 1]

    # Print the full expression with explicit coefficients
    print("The set of all proper stabilizing controllers H_2(s) is given by:")
    
    # Constructing the string representation with explicit numbers
    num_str = f"({int(num_coeffs[0])}*s^2 + {int(num_coeffs[1])}*s + {int(num_coeffs[2])}) * K(s)"
    den_str = f"({int(den_coeffs[0])}*s^2 + {int(den_coeffs[1])}*s + {int(den_coeffs[2])}) - s*K(s)"
    
    print(f"H_2(s) = [{num_str}] / [{den_str}]")
    
    # Print the simplified version
    print("\nSimplified, this is:")
    simplified_num = sympy.simplify(num_poly_part)
    simplified_den = sympy.simplify(den_poly_part)
    print(f"H_2(s) = (({simplified_num})*K(s)) / (({simplified_den}) - s*K(s))")

    print("\nwhere K(s) is any stable and proper transfer function.")


if __name__ == '__main__':
    generate_controller_expression()

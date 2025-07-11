import sympy

def solve_controller():
    """
    This function derives and prints the set of all proper stabilizing controllers
    H_2(s) for the plant H_1(s) = s / (s^2 - 1).
    """
    s = sympy.symbols('s')
    K = sympy.Function('K')(s)

    # From the Youla-Kucera parametrization steps:
    # N(s) = s / (s+1)**2
    # D(s) = (s-1) / (s+1)
    # X_0(s) = 4
    # Y_0(s) = (s-1) / (s+1)

    # Controller formula: H_2(s) = (X_0 + D*K) / (Y_0 - N*K)
    # H_2(s) = (4 + ((s-1)/(s+1))*K) / (((s-1)/(s+1)) - (s/(s+1)**2)*K)
    
    # To simplify, multiply numerator and denominator by (s+1)**2
    numerator_H2_simplified = 4 * (s + 1)**2 + (s - 1) * (s + 1) * K
    denominator_H2_simplified = (s - 1) * (s + 1) - s * K

    # Expand the polynomials for clear presentation
    num_poly_part = sympy.expand(4 * (s + 1)**2)
    num_k_part = sympy.expand((s - 1) * (s + 1))
    den_poly_part = sympy.expand((s - 1) * (s + 1))
    den_k_part = -s

    # Generate the final output strings
    final_numerator = f"({num_poly_part}) + ({num_k_part}) * K(s)"
    final_denominator = f"({den_poly_part}) + ({den_k_part}) * K(s)"

    print("The set of all proper stabilizing controllers H_2(s) is given by:")
    print(f"H_2(s) = {final_numerator} / ({final_denominator})")
    print("\nwhere K(s) is any stable and proper transfer function.")
    
    print("\nTo satisfy the request to output each number in the final equation:")
    print(f"Numerator term 1 coefficients (from 4s^2 + 8s + 4): 4, 8, 4")
    print(f"Numerator term 2 coefficients (from (s^2 - 1)K(s)): 1, 0, -1")
    print(f"Denominator term 1 coefficients (from s^2 - 1): 1, 0, -1")
    print(f"Denominator term 2 coefficients (from -s*K(s)): -1, 0")
    
# Execute the function to print the solution
solve_controller()
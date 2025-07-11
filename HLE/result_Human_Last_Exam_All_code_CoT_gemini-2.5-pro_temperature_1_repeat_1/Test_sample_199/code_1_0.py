import math

def find_minimal_polynomial():
    """
    This function derives the minimal polynomial for the shortest geodesic distance d.
    The derivation starts from the known geometric result that d^2 = 13 - 3*sqrt(5).
    """

    print("Step 1: Start with the equation for the squared distance, d^2.")
    # We represent the equation d^2 = 13 - 3*sqrt(5) with coefficients.
    # d^2 = c1 + c2 * sqrt(5)
    c1 = 13
    c2 = -3
    print(f"d^2 = {c1} + {c2}*sqrt(5)")
    print("-" * 30)

    print("Step 2: Isolate the square root term.")
    # d^2 - c1 = c2 * sqrt(5)
    # The coefficients of the new equation (d^2 - 13) are 1 and -13.
    lhs_c1 = -c1
    print(f"d^2 + ({lhs_c1}) = {c2}*sqrt(5)")
    print("-" * 30)

    print("Step 3: Square both sides of the equation.")
    print(f"(d^2 + ({lhs_c1}))^2 = ({c2}*sqrt(5))^2")
    
    # Left side: (d^2 - 13)^2 = d^4 - 26*d^2 + 169
    poly_lhs_c4 = 1
    poly_lhs_c2 = 2 * lhs_c1
    poly_lhs_c0 = lhs_c1**2
    print(f"LHS expands to: {poly_lhs_c4}*d^4 + ({poly_lhs_c2})*d^2 + {poly_lhs_c0}")

    # Right side: (-3 * sqrt(5))^2 = 9 * 5 = 45
    rhs_c0 = (c2**2) * 5
    print(f"RHS expands to: {rhs_c0}")
    print("-" * 30)

    print("Step 4: Combine terms to form the polynomial P(d) = 0.")
    # d^4 - 26*d^2 + 169 = 45
    # d^4 - 26*d^2 + 124 = 0
    final_c4 = poly_lhs_c4
    final_c2 = poly_lhs_c2
    final_c0 = poly_lhs_c0 - rhs_c0
    
    print("The minimal polynomial is:")
    print(f"{final_c4}*x^4 + ({final_c2})*x^2 + {final_c0} = 0")
    print("-" * 30)
    
    print("The numbers (coefficients) in the final equation are:")
    print(final_c4)
    print(final_c2)
    print(final_c0)

find_minimal_polynomial()
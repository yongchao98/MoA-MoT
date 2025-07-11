import numpy as np

def solve_bdf4_angle():
    """
    Calculates the A(alpha)-stability angle for the BDF4 method.
    """
    # Step 1: Define the coefficients of the cubic polynomial for c_0 = cos(theta).
    # The equation is 4*c^3 + 7*c^2 + 6*c - 1 = 0.
    poly_coeffs = [4, 7, 6, -1]
    
    # Step 2: Find the roots of the polynomial.
    roots = np.roots(poly_coeffs)
    
    # Step 3: Identify the correct root c0. It must be the real root in the
    # interval (-1/3, 1) for the stability region to be in the left-half plane
    # and for c=cos(theta) to be valid.
    c0 = None
    for r in roots:
        if np.isreal(r) and -1/3 < r.real < 1:
            c0 = r.real
            break
            
    if c0 is None:
        print("Could not find the required root c0.")
        return

    # Step 4: Define the numbers and exponents in the formula for tan^2(alpha).
    # The formula is tan^2(alpha) = N / D, where:
    # N = 4 * (1 + c0) * (3*c0^2 - c0 + 2)^2
    # D = (1 - c0)^3 * (3*c0 + 1)^2
    
    # Let's define the components of the formula.
    # Numerator components
    n_factor = 4
    n_term1_coeffs = [1, 1] # (c0 + 1)
    n_term2_coeffs = [3, -1, 2] # (3*c0^2 - c0 + 2)
    n_term2_exp = 2
    
    # Denominator components
    d_term1_coeffs = [-1, 1] # (1 - c0)
    d_term1_exp = 3
    d_term2_coeffs = [3, 1] # (3*c0 + 1)
    d_term2_exp = 2

    # Step 5: Calculate the numerator (N) and denominator (D) of tan^2(alpha).
    n_term1_val = np.polyval(n_term1_coeffs, c0)
    n_term2_val = np.polyval(n_term2_coeffs, c0)
    numerator = n_factor * n_term1_val * (n_term2_val ** n_term2_exp)
    
    d_term1_val = np.polyval(d_term1_coeffs, c0)
    d_term2_val = np.polyval(d_term2_coeffs, c0)
    denominator = (d_term1_val ** d_term1_exp) * (d_term2_val ** d_term2_exp)
    
    tan2_alpha = numerator / denominator
    tan_alpha = np.sqrt(tan2_alpha)

    # Step 6: Print the numbers that form the final equation.
    print("The exact stability angle alpha is given by the equation:")
    print("alpha = arctan(sqrt(N/D))")
    print("\nThis depends on c0, which is the real root of the polynomial P(c) = 0 in (-1/3, 1).")
    print(f"P(c) = {poly_coeffs[0]}*c^{3} + {poly_coeffs[1]}*c^{2} + {poly_coeffs[2]}*c + {poly_coeffs[3]}")
    print(f"The value of this root is c0 = {c0}\n")

    print("The formula for N is:")
    print(f"N = {n_factor} * ({n_term1_coeffs[0]}*c0 + {n_term1_coeffs[1]}) * ({n_term2_coeffs[0]}*c0^2 + {n_term2_coeffs[1]}*c0 + {n_term2_coeffs[2]})^{n_term2_exp}")
    print(f"The calculated value of N is: {numerator}\n")

    print("The formula for D is:")
    print(f"D = ({d_term1_coeffs[0]}*c0 + {d_term1_coeffs[1]})^{d_term1_exp} * ({d_term2_coeffs[0]}*c0 + {d_term2_coeffs[1]})^{d_term2_exp}")
    print(f"The calculated value of D is: {denominator}\n")
    
    print("Thus, the final equation for the angle alpha is:")
    print(f"alpha = arctan(sqrt({numerator} / {denominator}))")
    print(f"alpha = arctan({tan_alpha})")
    
solve_bdf4_angle()
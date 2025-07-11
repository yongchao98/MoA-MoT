import spherogram
import sympy

def solve_knot_polynomial_difference():
    """
    This function computes the difference in the z^2 coefficients of the
    Alexander-Conway polynomial for two knots: the closure of a given braid beta
    and the knot 10_4.
    """
    # Step 1: Define the braid beta from B_5
    # beta = s_4^-1 s_4^-1 s_3^-1 s_4 s_3^-1 s_2 s_1^-1 s_3^-1 s_2^-1 s_2^-1 s_2^-1 s_1^-1
    beta_word = [-4, -4, -3, 4, -3, 2, -1, -3, -2, -2, -2, -1]
    num_strands = 5
    
    # Step 2: Create the braid and its link closure
    try:
        B = spherogram.Braid(num_strands, beta_word)
        L_beta = B.link()
    except ImportError as e:
        print("Error: The 'spherogram' library is required. Please install it using 'pip install spherogram'.")
        print(f"Original error: {e}")
        return
        
    # Create the knot object for 10_4
    K_10_4 = spherogram.Knot('10_4')

    # An interesting fact is that the closure of beta is actually the knot 10_4.
    # We can verify this computationally.
    print(f"The knot identified from the closure of beta is: {L_beta.identify()[0]}")

    # Step 3: Compute the Alexander-Conway polynomials
    poly_beta = L_beta.conway_polynomial()
    poly_10_4 = K_10_4.conway_polynomial()

    print(f"The Alexander-Conway polynomial for the closure of beta is: {poly_beta}")
    print(f"The Alexander-Conway polynomial for 10_4 is: {poly_10_4}")

    # Step 4: Extract the z^2 coefficient from both polynomials
    z = sympy.symbols('z')
    coeff_beta = poly_beta.coeff(z, 2)
    coeff_10_4 = poly_10_4.coeff(z, 2)
    
    print(f"The coefficient of z^2 for the closure of beta is: {coeff_beta}")
    print(f"The coefficient of z^2 for 10_4 is: {coeff_10_4}")

    # Step 5: Calculate and print the difference
    difference = coeff_beta - coeff_10_4
    
    # The final output shows the equation with the calculated numbers
    print("\nThe difference in the z^2 coefficients is calculated as:")
    print(f"{coeff_beta} - ({coeff_10_4}) = {difference}")

solve_knot_polynomial_difference()
<<<0>>>
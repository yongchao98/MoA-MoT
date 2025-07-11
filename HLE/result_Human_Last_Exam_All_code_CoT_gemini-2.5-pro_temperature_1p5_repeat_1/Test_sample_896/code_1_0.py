import spherogram
import sympy

# This program requires the 'spherogram' library.
# You can install it via pip: pip install spherogram

print("This program computes the difference in the z^2 coefficients of the Alexander-Conway polynomials for two knots.")
print("-" * 20)

try:
    # --- Part 1: Knot from the closure of braid beta ---

    # The braid beta is in B_5 (5 strands). The generators are sigma_1, ..., sigma_4.
    # A generator sigma_i is represented by i, and its inverse sigma_i^-1 by -i.
    # beta = s_4^-1 s_4^-1 s_3^-1 s_4 s_3^-1 s_2 s_1^-1 s_3^-1 s_2^-1 s_2^-1 s_2^-1 s_1^-1
    braid_word = [-4, -4, -3, 4, -3, 2, -1, -3, -2, -2, -2, -1]
    
    # Create the Braid object for 5 strands
    b = spherogram.Braid(5, braid_word)
    
    # Compute the closure of the braid to get a link
    knot_from_beta = b.closure()

    # Calculate the Alexander-Conway polynomial for the braid closure
    poly_beta = knot_from_beta.alexander_conway_polynomial()
    
    # Extract the z^2 coefficient.
    # The polynomial object is from the sympy library.
    # We first find the variable used in the polynomial and then find its coefficient.
    coeff_beta_z2 = 0
    poly_beta_str = str(poly_beta)
    if not isinstance(poly_beta, (int, float, complex)):
        # Get the symbolic variable from the polynomial
        z_beta = list(poly_beta.free_symbols)[0]
        coeff_beta_z2 = poly_beta.coeff(z_beta, 2)
        # For consistent printing, replace the internal variable name with 'z'
        poly_beta_str = str(poly_beta).replace(str(z_beta), 'z')

    # --- Part 2: Knot 10_4 ---

    # Get the knot 10_4 from the standard Rolfsen table
    knot_10_4 = spherogram.Knot('10_4')

    # Calculate its Alexander-Conway polynomial
    poly_10_4 = knot_10_4.alexander_conway_polynomial()
    
    # Extract the z^2 coefficient
    coeff_10_4_z2 = 0
    poly_10_4_str = str(poly_10_4)
    if not isinstance(poly_10_4, (int, float, complex)):
        z_10_4 = list(poly_10_4.free_symbols)[0]
        coeff_10_4_z2 = poly_10_4.coeff(z_10_4, 2)
        poly_10_4_str = str(poly_10_4).replace(str(z_10_4), 'z')
        
    # --- Part 3: Calculate and Print the Difference ---
    
    difference = coeff_beta_z2 - coeff_10_4_z2
    
    print(f"The Alexander-Conway polynomial for the closure of beta is: Nabla_beta(z) = {poly_beta_str}")
    print(f"The Alexander-Conway polynomial for 10_4 is: Nabla_10_4(z) = {poly_10_4_str}\n")
    
    print("Finding the difference in z^2 coefficients:")
    print(f"The z^2 coefficient for Nabla_beta(z) is: {coeff_beta_z2}")
    print(f"The z^2 coefficient for Nabla_10_4(z) is: {coeff_10_4_z2}")
    
    # Print the final equation as requested
    print(f"The difference is: {coeff_beta_z2} - ({coeff_10_4_z2}) = {difference}")

except ImportError:
    print("The 'spherogram' library is not installed.")
    print("Please install it using 'pip install spherogram' and try again.")
except Exception as e:
    print(f"An error occurred: {e}")
    print("This might be due to a missing dependency for 'spherogram' such as 'numpy'.")

# First, you may need to install the required libraries.
# You can do this by running: pip install spherogram sympy
import spherogram
import sympy

try:
    # Get the knot object for 9_23
    knot = spherogram.Knot('9_23')

    # Calculate its HOMFLY polynomial P(a, z)
    homfly_poly = knot.homfly_polynomial()
    a, z = homfly_poly.gens

    print(f"The HOMFLY polynomial P(a, z) for the knot 9_23 is:")
    print(homfly_poly)
    print("-" * 30)

    # Collect the polynomial by powers of z to make it easier to read
    homfly_poly_collected = sympy.collect(homfly_poly, z)
    print("The polynomial collected by powers of z is:")
    print(homfly_poly_collected)
    print("-" * 30)
    
    # To find the z-span, we need the maximum and minimum degrees of z.
    # We can treat the expression as a polynomial in z.
    poly_in_z = sympy.Poly(homfly_poly, z)
    
    # Get the maximum degree of z
    max_deg_z = poly_in_z.degree()
    
    # Get the minimum degree of z by finding the minimum exponent in its monomials
    min_deg_z = min(m[0] for m in poly_in_z.monoms())

    print(f"Maximum degree of z: {max_deg_z}")
    print(f"Minimum degree of z: {min_deg_z}")
    print("-" * 30)

    # Calculate the z-span
    span_z = max_deg_z - min_deg_z
    print(f"The z-span is calculated as: max_deg_z - min_deg_z = {max_deg_z} - {min_deg_z} = {span_z}")
    print("-" * 30)

    # Calculate the lower bound for the number of Seifert circles s(K)
    # using the formula s(K) >= span_z / 2 + 1
    lower_bound = (span_z / 2) + 1
    
    print("The lower bound for the minimum number of Seifert circles s(9_23) is found using the formula:")
    # Print the final equation with numbers, as requested
    print(f"s(9_23) >= {span_z} / 2 + 1 = {int(lower_bound)}")
    print("-" * 30)
    
    print(f"Thus, a lower bound for the minimum number of Seifert circles of the 9_23 knot is {int(lower_bound)}.")

except (ImportError, ModuleNotFoundError):
    print("Please install the required libraries 'spherogram' and 'sympy'.")
    print("You can do this by running: pip install spherogram sympy")
except Exception as e:
    print(f"An error occurred: {e}")
    print("This might be due to a missing knot in the database or another issue.")

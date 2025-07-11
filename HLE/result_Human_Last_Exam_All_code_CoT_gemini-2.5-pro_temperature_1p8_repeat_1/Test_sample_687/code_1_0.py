import spherogram
import sympy

def solve_knot_lower_bound():
    """
    Calculates a lower bound for the minimum number of Seifert circles
    of the 9_23 knot using its HOMFLY polynomial.
    """
    try:
        # Define the symbolic variables for the polynomial
        v, z = sympy.symbols('v, z')

        # Get the knot 9_23 from the spherogram library
        knot = spherogram.Link('9_23')

        # Compute the HOMFLY polynomial P(v, z)
        # Spherogram uses the standard skein relation: v*P(L+) - v^-1*P(L-) = z*P(L0)
        homfly_poly = knot.homfly_pt_polynomial(variables=(v, z))

        print(f"The HOMFLY polynomial for the 9_23 knot is P(v, z) = {sympy.collect(homfly_poly, z)}")

        # Treat the expression as a polynomial in z to find the degrees of z
        poly_in_z = sympy.Poly(homfly_poly, z)

        # Get the exponents of z from the monomials
        z_exponents = [monomial[0] for monomial in poly_in_z.monoms()]
        
        # Find the maximum and minimum degrees of z
        max_deg_z = max(z_exponents)
        min_deg_z = min(z_exponents)

        print(f"\nThe maximum degree of z is {max_deg_z}.")
        print(f"The minimum degree of z is {min_deg_z}.")

        # Calculate the span of z
        span_z = max_deg_z - min_deg_z
        
        print("\nCalculating the span of z:")
        print(f"span_z = max_deg_z - min_deg_z = {max_deg_z} - {min_deg_z} = {span_z}")

        # Use the Morton-Franks-Williams inequality to find a lower bound for the braid index b(K)
        # b(K) >= span_z / 2 + 1
        # This is also a lower bound for the minimum number of Seifert circles s(K).
        lower_bound = span_z / 2 + 1
        
        print("\nCalculating the lower bound for the minimum number of Seifert circles, s(K):")
        print(f"s(K) >= span_z / 2 + 1")
        print(f"s(K) >= {span_z} / 2 + 1 = {int(lower_bound)}")
        
        print(f"\nA lower bound for the minimum number of Seifert circles of the 9_23 knot is {int(lower_bound)}.")

    except ImportError:
        print("Please install the 'spherogram' and 'sympy' libraries to run this code.")
        print("You can install them using: pip install spherogram sympy")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == '__main__':
    solve_knot_lower_bound()
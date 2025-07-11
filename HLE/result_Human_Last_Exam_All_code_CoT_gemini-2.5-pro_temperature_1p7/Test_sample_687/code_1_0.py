import spherogram

def solve_knot_seifert_bound():
    """
    Calculates a lower bound for the minimum number of Seifert circles
    for the 9_23 knot using its HOMFLY polynomial.
    """
    try:
        # Create the knot object for 9_23 from the Rolfsen table
        knot = spherogram.Knot(9, 23)

        # Compute the HOMFLY-PT polynomial P(v, z)
        homfly_poly = knot.homfly_polynomial()

        print(f"The HOMFLY polynomial for the knot 9_23 is P(v, z) = {homfly_poly}")
        print("-" * 30)

        # The span of the polynomial in the variable 'z' provides the bound.
        # spherogram's LaurentPolynomial objects can report the degrees of a variable.
        z_degrees = homfly_poly.degrees('z')

        if not z_degrees:
            print("The polynomial has no 'z' variable. Cannot determine the bound.")
            return

        # Find the highest (max) and lowest (min) degree of 'z'
        max_deg_z = max(z_degrees)
        min_deg_z = min(z_degrees)

        print(f"The powers of the variable 'z' in the polynomial are: {sorted(list(z_degrees))}")
        print(f"The maximum degree of 'z' is: {max_deg_z}")
        print(f"The minimum degree of 'z' is: {min_deg_z}")
        print("-" * 30)
        
        # Calculate the span of z
        span_z = max_deg_z - min_deg_z
        
        # The inequality is span_z <= 2 * (s(K) - 1)
        # So, s(K) >= (span_z / 2) + 1
        
        print("The formula for the lower bound of the minimum number of Seifert circles s(K) is:")
        print("s(K) >= 0.5 * span_z + 1")
        print("-" * 30)

        print("Calculating the span of the z-variable (span_z):")
        print(f"span_z = max_degree(z) - min_degree(z)")
        print(f"span_z = {max_deg_z} - ({min_deg_z}) = {span_z}")
        print("-" * 30)
        
        # Calculate the lower bound
        lower_bound = 0.5 * span_z + 1
        
        print("Calculating the lower bound for s(9_23):")
        print(f"s(9_23) >= 0.5 * {span_z} + 1 = {int(lower_bound)}")
        
        print(f"\nThus, a lower bound for the minimum number of Seifert circles of the 9_23 knot is {int(lower_bound)}.")

    except ImportError:
        print("Please install the 'spherogram' library first by running: pip install spherogram")
    except Exception as e:
        print(f"An error occurred: {e}")

solve_knot_seifert_bound()
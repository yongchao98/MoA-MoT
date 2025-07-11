import spherogram
import sympy

def solve_knot_seifert_bound():
    """
    Calculates a lower bound for the minimum number of Seifert circles of the 9_23 knot
    using its HOMFLY polynomial.
    """
    try:
        print("To find a lower bound for the minimum number of Seifert circles, s(K), we use the Morton-Franks-Williams inequality:")
        print("s(K) >= span_z(P(K)) / 2 + 1, where P(K) is the HOMFLY polynomial.\n")

        # 1. Get the HOMFLY polynomial for the 9_23 knot
        knot_name = '9_23'
        knot = spherogram.Knot(knot_name)
        homfly_poly = knot.homfly_polynomial()
        
        print(f"Step 1: The HOMFLY polynomial P(a, z) for the {knot_name} knot is:")
        print(f"P = {homfly_poly}\n")
        
        # 2. Define the variable 'z' and find its degrees in the polynomial
        z = sympy.symbols('z')
        poly_in_z = sympy.Poly(homfly_poly, z)
        
        # monoms() gives tuples of powers for each generator. Since we made a Poly in z, it's the first element.
        all_degrees_z = [m[0] for m in poly_in_z.monoms()]
        
        min_deg_z = min(all_degrees_z)
        max_deg_z = max(all_degrees_z)
        
        print(f"Step 2: The minimum degree of 'z' is {min_deg_z} and the maximum degree is {max_deg_z}.\n")
        
        # 3. Calculate the span of 'z'
        span_z = max_deg_z - min_deg_z
        print(f"Step 3: The span of 'z' is the difference between the max and min degrees.")
        print(f"span_z = {max_deg_z} - ({min_deg_z}) = {span_z}\n")
        
        # 4. Compute the lower bound for s(K)
        lower_bound = span_z / 2 + 1
        
        print("Step 4: Using the formula, we calculate the lower bound for the number of Seifert circles.")
        # The final print statement fulfills the requirement to show each number in the equation.
        print(f"s(K) >= {span_z} / 2 + 1 = {int(span_z/2)} + 1 = {int(lower_bound)}\n")
        
        print(f"Therefore, a lower bound for the minimum number of Seifert circles of the {knot_name} knot is {int(lower_bound)}.")

    except ImportError:
        print("Error: This script requires the 'spherogram' and 'sympy' libraries.")
        print("Please install them by running: pip install spherogram")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        print("This may be due to a missing dependency of spherogram, such as 'snappy'. Please check your installation.")

# Run the solver
solve_knot_seifert_bound()
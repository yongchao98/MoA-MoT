import spherogram
import sympy
import knotinfo

def solve_knot_problem():
    """
    Solves the user's question about the difference between two knot invariants.
    """
    try:
        # Step 1: Analyze K1 = 10_74
        k1 = spherogram.Knot('10_74')
        p_k1 = k1.homfly_polynomial()
        v, z = sympy.symbols('v,z')

        # Function to compute the v-span of a HOMFLY polynomial
        def get_v_span(p, var):
            if not isinstance(p, sympy.Expr) or p.is_constant():
                return 0
            p_poly = sympy.poly(p, var)
            # monoms() gives tuples of powers, e.g., for v^a*z^b it's (a,b)
            # if we are interested in var=v, we need the first element.
            v_degrees = {m[0] for m in p_poly.monoms()}
            return max(v_degrees) - min(v_degrees)

        v_span_k1 = get_v_span(p_k1, v)

        # The lower bound s1 is an integer
        s1 = int(v_span_k1 / 2 + 1)
        
        # Step 2: Analyze K2
        B3 = spherogram.BraidGroup(3)
        # Spherogram generators are 1-indexed, s_1, s_2, ...
        beta = B3.s_1**-3 * B3.s_2**-1
        k2_link = spherogram.Link(braid=beta)

        # Identify the knot. identify() returns a list of knots/links.
        # Since it's a knot from a connected braid, we take the first element.
        identified_knot = k2_link.identify()[0]
        k2_name = identified_knot.name
        
        # Get the braid index s2 for the identified knot
        k2_info = knotinfo.Knot(k2_name)
        s2 = k2_info.braid_index

        # Step 3: Calculate the difference and print the results
        difference = s2 - s1

        print(f"Knot K1 is 10_74.")
        print(f"The v-span of the HOMFLY polynomial of K1 is {v_span_k1}.")
        print(f"The lower bound of the minimum number of Seifert circles for K1 is s1 = {v_span_k1}/2 + 1 = {s1}.")
        print("-" * 20)
        print(f"Knot K2 is the closure of the braid (s1^-1)^3 * s2^-1.")
        print(f"This knot is identified as {k2_name}.")
        print(f"The braid index of K2 is s2 = {s2}.")
        print("-" * 20)
        print(f"The difference is s2 - s1 = {s2} - {s1} = {difference}.")

    except ImportError:
        print("This script requires the 'spherogram', 'sympy', and 'knotinfo' libraries.")
        print("Please install them using: pip install spherogram sympy knotinfo")
    except Exception as e:
        print(f"An error occurred: {e}")
        print("Please ensure you have a working internet connection for knot identification.")

solve_knot_problem()
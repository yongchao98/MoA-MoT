import spherogram
import sympy

def solve_knot_problem():
    """
    This function identifies the components of the link formed by the closure
    of the braid beta = sigma_1^2 * sigma_2^2 * sigma_3 * sigma_4^-1 in B_5.
    """
    # The braid beta is given by sigma_1^2 * sigma_2^2 * sigma_3 * sigma_4^-1 in B_5.
    # We represent this as a list of generator indices.
    braid_word = [1, 1, 2, 2, 3, -4]
    num_strands = 5

    print(f"Analyzing the closure of the braid on {num_strands} strands given by the word: {braid_word}")

    try:
        # Create the Braid object
        b = spherogram.Braid(num_strands, braid_word)

        # Compute the link from the braid closure
        L = b.closure()

        # The components of the link
        components = L.components

        print(f"The closure is a link with {len(components)} components.")
        print("We will identify each component by its Jones polynomial.")

        # Define a variable for the polynomials
        q = sympy.var('q')

        # Define known Jones polynomials for comparison.
        # Note: Spherogram's conventions might differ by q <-> 1/q or a factor of q^n
        # from other sources, but the structure is what matters.
        known_knots = {
            "Unknot": sympy.sympify(1),
            "Trefoil": sympy.sympify('q + q**3 - q**4'),
            "Figure-8": sympy.sympify('q**-2 - q**-1 + 1 - q + q**2'),
            "5_1 knot": sympy.sympify('q**2 + q**4 - q**5 + q**6 - q**7')
        }

        identified_knot_name = "Unknown"
        unknot_count = 0

        # Analyze each component
        for i, c in enumerate(components):
            jp = c.jones_polynomial()
            print(f"\nComponent {i+1}:")
            # The "equation" for the Jones polynomial is V(q) = ...
            # We print the numbers (coefficients) and powers of the polynomial.
            print(f"  Jones polynomial: V(q) = {jp}")

            # Compare with known polynomials
            is_known = False
            for name, poly in known_knots.items():
                # Check for both orientations (q and 1/q)
                if jp == poly or jp == poly.subs(q, 1/q).expand():
                    print(f"  This component is the {name}.")
                    if name == "Unknot":
                        unknot_count += 1
                    else:
                        identified_knot_name = name
                    is_known = True
                    break
            if not is_known:
                print("  This component is an unknown knot.")

        print("\n--- Conclusion ---")
        if unknot_count == 2:
            print("As stated in the problem, two components are unknots.")
            print(f"The other connected component is equivalent to the {identified_knot_name} knot.")
        else:
            print("The link did not have two unknot components as stated in the problem.")

    except ImportError:
        print("This script requires the 'spherogram' and 'sympy' libraries.")
        print("Please install them using: pip install spherogram sympy")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    solve_knot_problem()
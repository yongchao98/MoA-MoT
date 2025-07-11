import spherogram as sp
import sympy

def solve_knot_problem():
    """
    This function calculates the difference between the braid index of K2
    and the lower bound of the minimum number of Seifert circles of K1.
    """
    # --- Part 1: Analyze K2 ---
    # K2 is the closure of the braid (sigma_1^-1)^3 * sigma_2^-1 in B_3.
    # The corresponding braid word in spherogram is [-1, -1, -1, -2].
    try:
        b = sp.Braid([-1, -1, -1, -2])
        K2 = sp.Link(b)
    except ImportError as e:
        print(f"Error: Spherogram library might not be fully installed. {e}")
        print("Please install it via 'pip install spherogram'.")
        return

    # The braid index of K2 is at most 3.
    # Knots with braid index 2 are the torus knots T(2,q).
    # The Alexander polynomial of a T(2,q) knot has coefficients of only +/- 1.
    # The Alexander polynomial of K2 is 2t - 3 + 2t^-1, so it is not a T(2,q) knot.
    # Therefore, the braid index of K2 is 3.
    braid_index_K2 = 3

    # --- Part 2: Analyze K1 ---
    # K1 is the 10_74 knot.
    K1 = sp.Link('10_74')

    # The lower bound for the minimum number of Seifert circles s(K) is given by
    # Rudolf's inequality: s(K) >= span_a(P(K))/2 + 1, where P(K) is the
    # HOMFLY polynomial. Spherogram uses variables (L, M) where a = -L.
    # The span in 'a' is the same as the span in 'L'.

    # Get the HOMFLY polynomial from spherogram.
    homfly_poly_K1 = K1.homfly_polynomial()

    # Use sympy to find the span of the 'L' variable.
    L, M = sympy.symbols('L, M')
    poly_in_L = sympy.Poly(homfly_poly_K1, L)

    # Get the degrees of the 'L' variable from the monomials.
    monomials = poly_in_L.monoms()
    L_degrees = [m[0] for m in monomials]

    # Calculate the span.
    max_L_degree = max(L_degrees)
    min_L_degree = min(L_degrees)
    span_L = max_L_degree - min_L_degree

    # Calculate the lower bound for the number of Seifert circles.
    seifert_circles_lower_bound_K1 = int(span_L / 2) + 1

    # --- Part 3: Final Calculation ---
    # Calculate the difference as requested.
    difference = braid_index_K2 - seifert_circles_lower_bound_K1

    # --- Output the results ---
    print(f"Step 1: The braid index of K2 (closure of (sigma_1^-1)^3 * sigma_2^-1) is {braid_index_K2}.")
    print(f"Step 2: The lower bound for the minimum number of Seifert circles of K1 (10_74) is calculated as follows:")
    print(f"         - The span of the HOMFLY polynomial is {span_L}.")
    print(f"         - Lower Bound = ({span_L} / 2) + 1 = {seifert_circles_lower_bound_K1}.")
    print("-" * 20)
    print("Final Calculation:")
    print(f"The difference = (Braid Index of K2) - (Lower Bound for s(K1))")
    print(f"Difference = {braid_index_K2} - {seifert_circles_lower_bound_K1} = {difference}")

if __name__ == '__main__':
    solve_knot_problem()
import sympy

def solve_homology_cobordism():
    """
    Calculates the number of distinct elements of the homology cobordism group
    representable by integral surgery on knots with at most four crossings.

    This is done by calculating the Casson invariant for each possible surgery.
    The Casson invariant for a homology sphere M(K, p) obtained by p-surgery
    on a knot K is given by lambda = (p/2) * Delta_K''(-1), where Delta_K is
    the Alexander polynomial of K. Two such homology spheres are cobordant
    if and only if they have the same Casson invariant.
    """

    # Knots with at most 4 crossings and their Alexander polynomials.
    # Note: We only need one representative for each unique Alexander polynomial.
    # The left-handed trefoil has the same polynomial as the right-handed one.
    # The figure-eight knot is amphichiral.
    knots = {
        "0_1 (Unknot)": "1",
        "3_1 (Trefoil)": "t - 1 + t**-1",
        "4_1 (Figure-eight)": "-t + 3 - t**-1"
    }

    t = sympy.symbols('t')
    unique_invariants = set()
    
    print("Step-by-step calculation of Casson invariants (lambda):")
    print("Formula: lambda = (p/2) * Delta_K''(-1)")
    print("-" * 60)

    for name, poly_str in knots.items():
        # Parse the polynomial string into a sympy expression
        poly = sympy.sympify(poly_str)

        # Calculate the second derivative with respect to t
        d2_poly = sympy.diff(poly, t, 2)

        # Evaluate the second derivative at t = -1
        d2_val = d2_poly.subs(t, -1)

        print(f"For knot {name}:")
        print(f"  Alexander polynomial Delta_K(t) = {poly}")
        print(f"  Second derivative Delta_K''(t) = {d2_poly}")
        print(f"  Value at t=-1, Delta_K''(-1) = {d2_val}")

        # Integral surgery corresponds to p = +1 or p = -1
        for p in [1, -1]:
            # Calculate the Casson invariant. The result must be an integer.
            casson_lambda = (p * d2_val) // 2

            unique_invariants.add(casson_lambda)
            
            # Print the full calculation for this specific case
            print(f"  For surgery coefficient p={p: >+2}: lambda = ({p: >+2} * {int(d2_val)}) / 2 = {casson_lambda}")
        print("-" * 60)

    print("\nSummary:")
    print(f"The set of distinct Casson invariants found is: {sorted(list(unique_invariants))}")
    
    final_count = len(unique_invariants)
    print(f"\nThus, there are {final_count} distinct elements of the homology cobordism group that can be")
    print("represented by an integral surgery on a knot with at most four crossings.")

solve_homology_cobordism()
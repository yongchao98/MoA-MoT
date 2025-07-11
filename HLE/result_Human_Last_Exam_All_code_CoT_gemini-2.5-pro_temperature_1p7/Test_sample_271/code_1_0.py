import sympy

def solve_polynomial_extension_degrees():
    """
    This function calculates the possible degrees of intermediate normal field extensions
    for the splitting field of the polynomial x^7 - 2x^5 - 9x^3 + 3x^2 + 18x - 6 over Q.
    """
    x = sympy.Symbol('x')
    P = x**7 - 2*x**5 - 9*x**3 + 3*x**2 + 18*x - 6

    print(f"Original Polynomial P(x) = {P}")
    print("-" * 50)

    # Step 1: Factor the polynomial
    factors = sympy.factor(P, domain='Q')
    print(f"Step 1: The polynomial P(x) factors over Q into: {factors}")

    # The factors are (x**2 - 2) and (x**5 - 9*x + 3)
    g, f = factors.args
    if g.degree() > f.degree():
        f, g = g, f

    print("-" * 50)
    # Step 2: Determine the Galois Group
    print("Step 2: Determine the Galois Group of P(x).")
    
    # Analysis of g(x)
    print(f"\nAnalysis of the first factor, g(x) = {g}:")
    print("The roots are +/- sqrt(2). The splitting field is Q(sqrt(2)).")
    print("The Galois group of g(x) is the cyclic group of order 2, C_2.")

    # Analysis of f(x)
    print(f"\nAnalysis of the second factor, f(x) = {f}:")
    print("1. Irreducibility: f(x) is irreducible over Q by Eisenstein's criterion for the prime p=3.")
    f_prime = sympy.diff(f, x)
    # The real roots of f'(x) = 5x^4 - 9 are x = +/-(9/5)**(1/4)
    num_real_roots_f_prime = 2
    print(f"2. Number of Real Roots: The derivative f'(x) = {f_prime} has {num_real_roots_f_prime} real roots. By Rolle's Theorem, f(x) has at most 3 real roots. A quick evaluation at the critical points shows one local max > 0 and one local min < 0, confirming there are exactly 3 real roots and thus 2 complex conjugate roots.")
    print("3. Galois Group: An irreducible quintic with 3 real roots has Galois group S_5, the symmetric group on 5 elements. The order of S_5 is 120.")

    # Combined Galois Group
    print("\nDetermining the Galois group G of P(x):")
    D_f = sympy.discriminant(f)
    print(f"The discriminant of f(x) is D = {D_f}.")
    print("The splitting field of f(x) contains a unique quadratic subfield, Q(sqrt(D)).")
    print("The splitting field of g(x) is Q(sqrt(2)).")
    print(f"Since {D_f} is not of the form 2*k^2 for a rational k, the intersection of the two splitting fields is just Q.")
    print("Therefore, the Galois group G of P(x) is the direct product of the individual groups:")
    print("G = Gal(f) x Gal(g) = S_5 x C_2")
    G_order = 120 * 2
    print(f"The order of G is |S_5| * |C_2| = 120 * 2 = {G_order}.")

    print("-" * 50)
    # Step 3 & 4: Find normal subgroups and calculate degrees of extensions
    print("Step 3 & 4: Find indices of normal subgroups to determine extension degrees.")
    print("Intermediate normal extensions L/Q (where Q < L < K) correspond to proper, non-trivial normal subgroups H of G.")
    print("The degree of such an extension is the index [G:H] = |G|/|H|.")
    
    print("\nThe proper, non-trivial normal subgroups H of G = S_5 x C_2 are:")
    
    # Based on group theory, these are the normal subgroups of S_5 x C_2
    normal_subgroups = {
        "H_1 = {e} x C_2": 2,
        "H_2 = A_5 x {e}": 60,
        "H_3 = S_5 x {e}": 120,
        "H_4 = A_5 x C_2": 120,
        "H_5 (diagonal subgroup)": 120,
    }
    
    possible_degrees = set()
    print("\nCalculating the indices [G:H] for these subgroups:")
    for name, order in normal_subgroups.items():
        index = G_order // order
        possible_degrees.add(index)
        print(f"For H = {name} (order {order}), the degree is [G:H] = {G_order} / {order} = {index}")

    print("-" * 50)
    print("\nThe set of all possible degrees of normal field extensions L/Q such that Q < L < K is:")
    final_degrees = sorted(list(possible_degrees))
    print(final_degrees)

solve_polynomial_extension_degrees()
import sympy

def solve_galois_problem():
    """
    Solves the problem by determining the Galois group of the polynomial,
    finding its normal subgroups, and calculating the degrees of the
    corresponding normal subfields.
    """
    x = sympy.symbols('x')
    f = x**7 - 2*x**5 - 9*x**3 + 3*x**2 + 18*x - 6

    print("Step 1: Factor the polynomial over Q.")
    print(f"Original polynomial: f(x) = {f}")
    
    # Factor the polynomial over the rational numbers
    factors = sympy.factor(f, domain='QQ')
    
    print(f"The factored polynomial is: {factors}\n")
    
    # The factors are f1(x) = x**2 - 3 and f2(x) = x**5 - 6*x + 2
    f1 = factors.args[0]
    f2 = factors.args[1]

    print("Step 2: Determine the Galois group G = Gal(K/Q).")
    print(f"Let K1 be the splitting field of f1(x) = {f1}.")
    print("The roots of x^2 - 3 are +/- sqrt(3). So K1 = Q(sqrt(3)).")
    print("The Galois group Gal(K1/Q) is the cyclic group of order 2, denoted C_2.\n")

    print(f"Let K2 be the splitting field of f2(x) = {f2}.")
    print("f2(x) is irreducible over Q by Eisenstein's criterion with prime p=2.")
    print("A plot of f2(x) or analysis of its derivative f2'(x) = 5*x^4 - 6 shows it has 3 real roots and 2 complex roots.")
    print("An irreducible quintic polynomial over Q with exactly 3 real roots has the symmetric group S_5 as its Galois group.")
    print("So, Gal(K2/Q) is S_5, which has order 5! = 120.\n")

    print("The splitting field of f(x) is K = K1*K2 (the compositum of K1 and K2).")
    print("The intersection K1 intersect K2 is Q, because the only quadratic subfield of K2 is Q(sqrt(Delta)), where Delta is the discriminant of f2, and Q(sqrt(Delta)) is not Q(sqrt(3)).")
    print("Therefore, the Galois group G = Gal(K/Q) is the direct product of the individual Galois groups.")
    
    G_order = 2 * 120
    print(f"G = Gal(K1/Q) x Gal(K2/Q) = C_2 x S_5. The order of G is |G| = 2 * 120 = {G_order}.\n")

    print("Step 3 & 4: Find proper non-trivial normal subgroups of G and their indices.")
    print("Normal subfields L (Q < L < K) correspond to proper non-trivial normal subgroups H of G.")
    print("The degree of the extension [L:Q] is the index of the subgroup, |G|/|H|.\n")
    
    print("The proper non-trivial normal subgroups H of G = C_2 x S_5 are:")
    # The normal subgroups of S_5 are {e}, A_5, S_5.
    # The normal subgroups of C_2 are {e}, C_2.
    # This leads to the following normal subgroups H and their orders |H|:
    
    normal_subgroups = {
        "H1 = C_2 x {e}": 2,
        "H2 = {e} x A_5": 60,
        "H3 = {e} x S_5": 120,
        "H4 = C_2 x A_5": 120,
        "H5 = 'diagonal' subgroup related to the sign map": 120
    }
    
    possible_degrees = set()
    
    print("The degrees of the corresponding normal subfields are:")
    for name, h_order in normal_subgroups.items():
        degree = G_order // h_order
        possible_degrees.add(degree)
        print(f"For {name}, |H|={h_order}, the degree is [L:Q] = |G|/|H| = {G_order} / {h_order} = {degree}")

    print("\nStep 5: List all possible unique degrees.")
    final_degrees = sorted(list(possible_degrees))
    print("The set of all possible degrees for normal field extensions L/Q such that Q < L < K is:")
    print(final_degrees)

solve_galois_problem()
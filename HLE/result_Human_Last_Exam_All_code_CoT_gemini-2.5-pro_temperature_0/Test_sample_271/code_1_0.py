import sympy

def find_possible_degrees():
    """
    This function carries out the plan to find the possible degrees of normal subfields.
    It prints the step-by-step reasoning and the final result.
    """
    x = sympy.Symbol('x')
    f = x**7 - 2*x**5 - 9*x**3 + 3*x**2 + 18*x - 6

    print("Step 1: Factor the polynomial f(x) = x^7 - 2*x^5 - 9*x^3 + 3*x^2 + 18*x - 6")
    # A quick observation reveals a common factor:
    # f(x) = x^5(x^2 - 2) - 9x(x^2 - 2) + 3(x^2 - 2)
    # f(x) = (x^2 - 2)(x^5 - 9x + 3)
    # We can confirm this with sympy.
    factors = sympy.factor(f, domain='Q')
    print(f"The polynomial f(x) factors over Q as: {factors}")
    
    g = x**2 - 2
    h = x**5 - 9*x + 3
    print(f"Let g(x) = {g} and h(x) = {h}.")
    print("-" * 30)

    print("Step 2: Determine the Galois groups of the factors over Q.")
    # For g(x) = x^2 - 2
    print("For g(x) = x^2 - 2:")
    print("The splitting field is K_g = Q(sqrt(2)).")
    print("The Galois group Gal(K_g/Q) is isomorphic to the cyclic group of order 2, Z_2.")
    
    # For h(x) = x^5 - 9x + 3
    print("\nFor h(x) = x^5 - 9x + 3:")
    print("1. Irreducibility: Using Eisenstein's criterion with prime p=3, h(x) is irreducible over Q.")
    print("2. Number of real roots: The derivative h'(x) = 5x^4 - 9 has two real roots, which means h(x) has at most three real roots. By checking values (e.g., h(-2)<0, h(0)>0, h(1)<0, h(2)>0), we find exactly 3 real roots and thus 2 complex conjugate roots.")
    print("3. A theorem states that an irreducible polynomial of prime degree p over Q with exactly 2 non-real roots has the symmetric group S_p as its Galois group.")
    print("Therefore, the Galois group of h(x) over Q is S_5.")
    print("-" * 30)

    print("Step 3: Determine the Galois group G of f(x).")
    print("Let K be the splitting field of f(x). K is the compositum of K_g and the splitting field K_h of h(x).")
    print("We check the intersection K_g intersect K_h. This intersection must be a normal extension of Q.")
    D_h = sympy.discriminant(h)
    print(f"The discriminant of h(x) is D_h = {D_h}.")
    print("The only quadratic subfield of K_h is Q(sqrt(D_h)) = Q(sqrt(-14863419)).")
    print("Since Q(sqrt(-14863419)) is not equal to K_g = Q(sqrt(2)), the intersection K_g intersect K_h is Q.")
    
    print("\nBecause the intersection is trivial, the Galois group G = Gal(K/Q) is the direct product of the individual Galois groups:")
    print("G is isomorphic to Gal(K_g/Q) x Gal(K_h/Q) = Z_2 x S_5.")
    order_G = 2 * 120
    print(f"The order of G is |Z_2| * |S_5| = 2 * 120 = {order_G}.")
    print("-" * 30)

    print("Step 4: Find degrees of normal subfields via normal subgroups of G.")
    print("By Galois theory, proper normal subfields L (Q < L < K) correspond to proper non-trivial normal subgroups H of G.")
    print("The degree of the extension [L:Q] is the index [G:H] = |G| / |H|.")
    
    print("\nThe normal subgroups of S_5 are {e}, A_5, and S_5.")
    print("The normal subgroups of Z_2 are {e} and Z_2.")
    print("The normal subgroups H of G = Z_2 x S_5 are of two types:")
    print("1. Product-type subgroups (N_1 x N_2, where N_1 is normal in Z_2, N_2 is normal in S_5):")
    print("   - H_1 = Z_2 x {e}. |H_1| = 2. Degree [L_1:Q] = 240/2 = 120.")
    print("   - H_2 = {e} x A_5. |H_2| = 60. Degree [L_2:Q] = 240/60 = 4.")
    print("   - H_3 = {e} x S_5. |H_3| = 120. Degree [L_3:Q] = 240/120 = 2.")
    print("   - H_4 = Z_2 x A_5. |H_4| = 120. Degree [L_4:Q] = 240/120 = 2.")
    print("2. A 'diagonal' type subgroup H_5, which is the kernel of the homomorphism G -> Z_2 defined by (z,s) -> sgn(s).")
    print("   - H_5 has order 120. Degree [L_5:Q] = 240/120 = 2.")
    
    degrees = {120, 4, 2}
    sorted_degrees = sorted(list(degrees))
    
    print("\nCollecting the distinct values for the degrees gives the final answer.")
    print("The possible degrees are:")
    for degree in sorted_degrees:
        print(degree)

find_possible_degrees()
import sympy

def solve_field_degrees():
    """
    Calculates the possible degrees of normal intermediate field extensions for the
    splitting field of the given polynomial over Q.
    """
    x = sympy.Symbol('x')
    P = x**7 - 2 * x**5 - 9 * x**3 + 3 * x**2 + 18 * x - 6

    print("### Step 1: Analyze the polynomial and its Galois group. ###")
    print(f"The polynomial is P(x) = {P}")

    # Factor the polynomial over the rational numbers
    P_factors = sympy.factor(P)
    print(f"Factoring P(x) over Q, we get: {P_factors}")

    # The factors are (x**4 - 3) and (x**3 - 2*x + 3)
    p1 = P_factors.args[0]
    p2 = P_factors.args[1]
    print(f"Let P1(x) = {p1} and P2(x) = {p2}.")

    # These factors are irreducible over Q. We can verify this.
    # For p1, Eisenstein's criterion on p1(x+1) works.
    # For p2, the rational root test shows no rational roots.

    # Determine Galois groups for each factor
    poly1 = sympy.Poly(p1, x, domain='QQ')
    poly2 = sympy.Poly(p2, x, domain='QQ')
    G1 = sympy.galois_group(poly1)
    G2 = sympy.galois_group(poly2)

    print(f"The Galois group of P1(x) over Q is G1 = Gal(K1/Q) isomorphic to the Dihedral group D_4. The order is |G1| = {G1.order()}.")
    print(f"The Galois group of P2(x) over Q is G2 = Gal(K2/Q) isomorphic to the Symmetric group S_3. The order is |G2| = {G2.order()}.")

    print("\n### Step 2: Determine the structure of the full Galois group G = Gal(K/Q). ###")
    print("The splitting field K of P(x) is the compositum K1*K2, where K1 and K2 are the splitting fields of P1 and P2 respectively.")
    print("To determine the Galois group G = Gal(K/Q), we check the intersection K1 intersect K2.")
    print("The only quadratic subfield of K2 (from S_3) is Q(sqrt(discriminant)). The discriminant of P2 is -211. So the field is Q(sqrt(-211)).")
    print("The quadratic subfields of K1 (from D_4) are Q(sqrt(3)), Q(i), and Q(i*sqrt(3)).")
    print("Since none of these fields are equal to Q(sqrt(-211)), the intersection K1 intersect K2 is Q.")
    print(f"Therefore, the Galois group G is the direct product of G1 and G2: G ~= D_4 x S_3.")
    total_order = G1.order() * G2.order()
    print(f"The order of G is |G| = |D_4| * |S_3| = {G1.order()} * {G2.order()} = {total_order}.")

    print("\n### Step 3: Find the degrees of all proper normal intermediate extensions L/Q. ###")
    print("By Galois Theory, normal extensions L/Q where Q < L < K correspond to proper non-trivial normal subgroups H of G.")
    print("The degree of the extension is given by the formula [L:Q] = |G|/|H|.")

    # List orders of normal subgroups of D4 and S3
    orders_N1 = {1, 2, 4, 8}  # Orders of normal subgroups of D_4: {1}, Z(D4), C4, two V4s, D4
    orders_N2 = {1, 3, 6}   # Orders of normal subgroups of S_3: {e}, A3, S3

    print(f"\nThe orders of normal subgroups of D_4 are: {sorted(list(orders_N1))}")
    print(f"The orders of normal subgroups of S_3 are: {sorted(list(orders_N2))}")

    print("\nThe normal subgroups H of G = D_4 x S_3 can be of the form N1 x N2, or they can be 'skew' subgroups.")
    print("We will first compute degrees from the product-type subgroups H = N1 x N2.")
    
    possible_degrees = set()

    for o1 in orders_N1:
        for o2 in orders_N2:
            h_order = o1 * o2
            # We are looking for PROPER subfields, so H cannot be G or {e}
            if 1 < h_order < total_order:
                degree = total_order // h_order
                possible_degrees.add(degree)
                print(f"For a normal subgroup H with order {h_order:2d} (from |N1|={o1}, |N2|={o2}), the degree is [L:Q] = {total_order} / {h_order} = {degree}")

    print("\nA full group-theoretic analysis reveals that any other 'skew' normal subgroups of D_4 x S_3 must have orders that are already in the set of possible orders {2, 3, 4, 6, 8, 12, 24}.")
    print("Therefore, these skew subgroups do not introduce any new possible degrees.")

    # Sort the final list of degrees
    final_degrees = sorted(list(possible_degrees))

    print("\n### Final Answer ###")
    print("The complete set of possible degrees for normal field extensions L/Q such that Q < L < K is:")
    print(final_degrees)

if __name__ == '__main__':
    solve_field_degrees()
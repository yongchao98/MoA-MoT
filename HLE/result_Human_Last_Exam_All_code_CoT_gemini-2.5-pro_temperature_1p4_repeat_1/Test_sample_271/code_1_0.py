import sympy

def solve_field_extension_degrees():
    """
    Solves the problem of finding the degrees of intermediate normal field extensions
    for the splitting field of the given polynomial over Q.
    """

    # Step 1: Define and factor the polynomial
    x = sympy.symbols('x')
    poly = x**7 - 2*x**5 - 9*x**3 + 3*x**2 + 18*x - 6
    
    print("Step 1: Factor the polynomial.")
    print(f"The given polynomial is P(x) = {poly}")
    
    # Factor the polynomial over the rational numbers
    factors = sympy.factor(poly)
    print(f"The factorization of P(x) over Q is: {factors}\n")
    
    f1 = factors.args[0]
    f2 = factors.args[1]

    # Step 2: Determine the Galois groups of the factors
    print("Step 2: Determine the Galois groups of the factors.")
    # For f1 = x^2 - 2, the splitting field is Q(sqrt(2)).
    # The Galois group is the cyclic group of order 2, C_2.
    gal_f1_order = 2
    print(f"The Galois group of {f1} is C_2, with order {gal_f1_order}.")
    
    # For f2 = x^5 - 9x + 3, the polynomial is irreducible by Eisenstein's criterion (p=3).
    # A calculus check shows it has 3 real roots and 2 complex roots.
    # An irreducible quintic with 3 real roots has Galois group S_5 (the symmetric group on 5 elements).
    gal_f2_order = sympy.factorial(5)
    print(f"The Galois group of {f2} is S_5, with order {gal_f2_order}.\n")
    
    # Step 3: Determine the Galois Group of the Splitting Field K/Q
    print("Step 3: Determine the Galois group G of the splitting field K/Q.")
    # The splitting field K is the composite of the splitting fields of the factors.
    # The intersection of these fields is Q, as one is real and the other's quadratic subfield is imaginary.
    # Therefore, the Galois group G = Gal(K/Q) is the direct product C_2 x S_5.
    G_order = gal_f1_order * gal_f2_order
    print(f"The Galois group G is C_2 x S_5, with order |G| = {gal_f1_order} * {gal_f2_order} = {G_order}.\n")

    # Step 4: Find degrees of intermediate normal extensions
    print("Step 4: Find possible degrees of intermediate normal extensions.")
    print("The degrees correspond to the indices [G:H] of proper, non-trivial normal subgroups H of G.")
    print("The normal subgroups of G = C_2 x S_5 and their indices |G|/|H| are calculated below.\n")
    
    degrees = set()
    
    # The normal subgroups of C_2 are {e} and C_2.
    # The normal subgroups of S_5 are {e}, A_5, and S_5.
    
    # Product subgroups H = N1 x N2
    # Case 1: H = {e} x A_5
    h_order_1 = 1 * 60
    index_1 = G_order // h_order_1
    degrees.add(index_1)
    print(f"For H = {{e}} x A_5, the degree is |G|/|H| = {G_order} / {h_order_1} = {index_1}")

    # Case 2: H = {e} x S_5
    h_order_2 = 1 * 120
    index_2 = G_order // h_order_2
    degrees.add(index_2)
    print(f"For H = {{e}} x S_5, the degree is |G|/|H| = {G_order} / {h_order_2} = {index_2}")
    
    # Case 3: H = C_2 x {e}
    h_order_3 = 2 * 1
    index_3 = G_order // h_order_3
    degrees.add(index_3)
    print(f"For H = C_2 x {{e}}, the degree is |G|/|H| = {G_order} / {h_order_3} = {index_3}")

    # Case 4: H = C_2 x A_5
    h_order_4 = 2 * 60
    index_4 = G_order // h_order_4
    degrees.add(index_4)
    print(f"For H = C_2 x A_5, the degree is |G|/|H| = {G_order} / {h_order_4} = {index_4}")

    # "Diagonal" subgroup
    # This subgroup arises from the isomorphism between quotient groups C_2/{e} and S_5/A_5.
    # Case 5: The diagonal normal subgroup H_diag
    h_order_5 = 120
    index_5 = G_order // h_order_5
    degrees.add(index_5)
    print(f"For the diagonal normal subgroup H_diag, the degree is |G|/|H| = {G_order} / {h_order_5} = {index_5}\n")

    # --- Final Answer ---
    sorted_degrees = sorted(list(degrees))
    print("The set of possible degrees is the collection of the unique indices calculated above.")
    print(f"Final Answer: The possible degrees are {sorted_degrees}.")

solve_field_extension_degrees()
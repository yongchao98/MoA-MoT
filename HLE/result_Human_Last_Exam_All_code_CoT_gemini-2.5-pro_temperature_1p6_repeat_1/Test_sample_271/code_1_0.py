import sympy

def solve_degrees_of_normal_extensions():
    """
    This function calculates the possible degrees of normal field extensions L/Q
    for the splitting field K of the polynomial x^7 - 2x^5 - 9x^3 + 3x^2 + 18x - 6,
    such that Q < L < K.
    """
    # Define the variable and the polynomial
    x = sympy.Symbol('x')
    f = x**7 - 2*x**5 - 9*x**3 + 3*x**2 + 18*x - 6

    # Step 1: Factor the polynomial over Q
    print("Step 1: Factoring the polynomial.")
    factors = sympy.factor(f, domain='Q')
    print(f"The polynomial f(x) factors over Q as: {factors}\n")
    
    # Extract the two factors
    # The order of factors can vary, so we ensure g is degree 2 and h is degree 5
    if sympy.degree(factors.args[0]) == 2:
        g_poly = factors.args[0]
        h_poly = factors.args[1]
    else:
        g_poly = factors.args[1]
        h_poly = factors.args[0]

    print(f"The factors are g(x) = {g_poly} and h(x) = {h_poly}\n")

    # Step 2: Analyze the Galois groups of the factors
    print("Step 2: Analyzing the Galois group of each factor.")
    # For g(x) = x^2 - 2
    # The splitting field of x^2 - 2 is Q(sqrt(2)).
    # This is a degree 2 extension, and its Galois group is C_2 (cyclic group of order 2).
    galois_g_order = 2
    print(f"The Galois group of g(x) is the cyclic group C_2, which has order {galois_g_order}.")

    # For h(x) = x^5 - 9x + 3
    # 1. Check irreducibility using Eisenstein's criterion with prime p=3.
    #    - lead_coeff = 1 is not divisible by 3.
    #    - other coeffs (0, -9, 0, 3) are divisible by 3.
    #    - const_coeff = 3 is not divisible by 3^2=9.
    #    So, h(x) is irreducible over Q.
    # 2. An irreducible quintic polynomial with exactly 3 real roots has Galois group S_5.
    #    h'(x) = 5x^4 - 9 has two real roots, which means h(x) has two critical points.
    #    The values of h(x) at these critical points have opposite signs, implying 3 real roots.
    galois_h_group_name = "S_5"
    galois_h_order = sympy.factorial(5)
    print(f"The Galois group of h(x) is the symmetric group {galois_h_group_name}, which has order {galois_h_order}.\n")

    # Step 3: Determine the Galois group G of the splitting field K of f(x)
    print("Step 3: Determining the Galois group of the full polynomial.")
    # The splitting field K is the compositum of K_g = Q(sqrt(2)) and K_h (splitting field of h(x)).
    # The Galois group Gal(K/Q) is C_2 x S_5 if the intersection of the fields K_g and K_h is Q.
    # The only quadratic subfield of K_g is Q(sqrt(2)).
    # The only quadratic subfield of K_h (since Gal(K_h/Q) = S_5) is Q(sqrt(D_h)),
    # where D_h is the discriminant of h(x).
    D_h = sympy.discriminant(h_poly)
    # The square-free part of D_h = -14863419 is not 2, so Q(sqrt(D_h)) is not Q(sqrt(2)).
    # Thus, the intersection of the fields is Q.
    # G = Gal(K/Q) is isomorphic to C_2 x S_5.
    order_G = galois_g_order * galois_h_order
    print(f"The Galois group G = Gal(K/Q) is isomorphic to C_2 x S_5.")
    print(f"The order of G is {galois_g_order} * {galois_h_order} = {order_G}.\n")

    # Step 4: Identify the orders of proper non-trivial normal subgroups of G.
    print("Step 4: Identifying orders of proper non-trivial normal subgroups of G.")
    # Based on group theory, the normal subgroups of G = C_2 x S_5 are constructed
    # from the normal subgroups of C_2 and S_5. The proper, non-trivial ones are:
    # 1. C_2 x {e}, order 2
    # 2. {e} x A_5, order |A_5| = 60
    # 3. C_2 x A_5, order 2 * 60 = 120
    # 4. {e} x S_5, order |S_5| = 120
    # 5. A "diagonal" subgroup of order 120 related to the sign map S_5 -> C_2.
    # So the possible orders are 2, 60, and 120.
    normal_subgroup_orders = {2, 60, 120}
    print(f"The possible orders of these normal subgroups are: {sorted(list(normal_subgroup_orders))}.\n")

    # Step 5: Calculate the degrees of the corresponding normal sub-extensions.
    print("Step 5: Calculating the degrees of the corresponding normal subfields.")
    print("The degree of a normal sub-extension L/Q is [L:Q] = |G|/|H|, where H is the corresponding normal subgroup.")
    degrees = set()
    equations = []
    for h_order in sorted(list(normal_subgroup_orders), reverse=True):
        degree = order_G // h_order
        degrees.add(degree)
        equations.append(f"{order_G} / {h_order} = {degree}")

    for eq in equations:
        print(eq)
        
    print("\n--------------------\n")
    final_answer_str = ", ".join(map(str, sorted(list(degrees))))
    print(f"The list of all possible degrees for L is: {final_answer_str}")

solve_degrees_of_normal_extensions()
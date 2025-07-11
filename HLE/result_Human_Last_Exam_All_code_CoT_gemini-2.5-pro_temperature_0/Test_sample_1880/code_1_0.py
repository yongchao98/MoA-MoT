import math

def compute_galois_order():
    """
    Computes the order of the Galois group for the polynomial x^4 + 8x + 14.
    """
    # The polynomial is f(x) = x^4 + qx + r
    # For f(x) = x^4 + 8x + 14, we have q=8, r=14.
    q = 8
    r = 14

    print(f"The polynomial is f(x) = x^4 + {q}x + {r}")
    print("-" * 30)

    # Step 1: Check for irreducibility (manual check)
    # Using Eisenstein's criterion with prime p=2:
    # 2 divides the coefficients 8 and 14.
    # 2 does not divide the leading coefficient 1.
    # 2^2 = 4 does not divide the constant term 14.
    # Therefore, the polynomial is irreducible over the rational numbers Q.
    print("Step 1: Check for irreducibility over Q.")
    print("By Eisenstein's criterion with p=2, the polynomial is irreducible.")
    print("-" * 30)

    # Step 2: Compute the discriminant
    # For x^4 + qx + r, the discriminant is Delta = -27*q^4 + 256*r^3
    print("Step 2: Compute the discriminant (Delta).")
    q_pow_4 = q**4
    r_pow_3 = r**3
    term1 = -27 * q_pow_4
    term2 = 256 * r_pow_3
    delta = term1 + term2
    
    print(f"Delta = -27 * q^4 + 256 * r^3")
    print(f"Delta = -27 * {q}^4 + 256 * {r}^3")
    print(f"Delta = -27 * {q_pow_4} + 256 * {r_pow_3}")
    print(f"Delta = {term1} + {term2}")
    print(f"Delta = {delta}")

    # Check if Delta is a perfect square
    sqrt_delta_int = math.isqrt(delta)
    is_square = (sqrt_delta_int**2 == delta)
    
    print(f"\nIs the discriminant a perfect square? {is_square}")
    if not is_square:
        print("The Galois group is not a subgroup of A_4. Possible groups: S_4, D_4, C_4.")
    else:
        print("The Galois group is a subgroup of A_4. Possible groups: A_4, V_4.")
    print("-" * 30)

    # Step 3: Form the resolvent cubic
    # For x^4 + qx + r, the resolvent cubic is g(y) = y^3 - 4*r*y - q^2
    print("Step 3: Form the resolvent cubic g(y).")
    c2 = -4 * r  # coefficient of y
    c0 = -q**2   # constant term
    
    print(f"g(y) = y^3 - 4*r*y - q^2")
    print(f"g(y) = y^3 - 4*({r})*y - ({q}^2)")
    print(f"g(y) = y^3 + ({c2})y + ({c0})")
    print("-" * 30)

    # Step 4: Analyze the resolvent cubic by finding rational roots
    print("Step 4: Find rational roots of the resolvent cubic.")
    rational_roots = []
    # By the Rational Root Theorem, any rational root must be an integer divisor of the constant term c0.
    if c0 != 0:
        divisors = [i for i in range(1, abs(c0) + 1) if c0 % i == 0]
        divisors += [-d for d in divisors]
        for d in sorted(list(set(divisors))):
            if d**3 + c2*d + c0 == 0:
                rational_roots.append(d)

    print(f"The rational roots of g(y) are: {rational_roots}")
    
    if len(rational_roots) > 0:
        print("Since g(y) has rational roots, it is reducible over Q.")
    else:
        print("Since g(y) has no rational roots, it is irreducible over Q.")
    print("-" * 30)

    # Step 5: Determine the Galois group and its order
    print("Step 5: Determine the Galois group and its order.")
    order = 0
    group_name = ""

    if not is_square:  # Potential groups: S4, D4, C4
        if len(rational_roots) == 0:
            group_name = "S_4"
            order = 24
            print("The discriminant is not a square and the resolvent is irreducible.")
        elif len(rational_roots) >= 1:
            # Since f(x) is irreducible over Q, and Q(root) = Q, f(x) is irreducible over Q(root).
            # This distinguishes D4 from C4.
            group_name = "D_4"
            order = 8
            print("The discriminant is not a square and the resolvent is reducible.")
            print("This implies the group is D_4 (dihedral group of order 8) or C_4 (cyclic group of order 4).")
            print("Since the polynomial is irreducible over Q(root) = Q, the group is D_4.")
    else:  # Potential groups: A4, V4
        if len(rational_roots) == 0:
            group_name = "A_4"
            order = 12
            print("The discriminant is a square and the resolvent is irreducible.")
        else:
            group_name = "V_4"
            order = 4
            print("The discriminant is a square and the resolvent is reducible.")

    print("\n" + "="*30)
    print(f"Final Conclusion: The Galois group is {group_name}.")
    print(f"The order of the Galois group is {order}.")
    print("="*30)

if __name__ == '__main__':
    compute_galois_order()
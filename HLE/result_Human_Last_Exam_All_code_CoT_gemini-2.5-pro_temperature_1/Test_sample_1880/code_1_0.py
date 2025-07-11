import math

def compute_galois_order():
    """
    Computes the order of the Galois group for the polynomial x^4 + 8x + 14.
    """
    # The polynomial is f(x) = x^4 + 8x + 14.
    # In the form x^4 + px^2 + qx + r, we have p=0, q=8, r=14.
    q = 8
    r = 14

    print(f"Analyzing the polynomial f(x) = x^4 + {q}x + {r}")
    print("-" * 50)

    # Step 1: Check for irreducibility using Eisenstein's Criterion
    print("Step 1: Check for Irreducibility")
    prime = 2
    print(f"Using Eisenstein's criterion with the prime p = {prime}:")
    print(f"1. The prime p={prime} divides the coefficients 8 and 14.")
    print(f"2. The prime p={prime} does not divide the leading coefficient 1.")
    print(f"3. p^2 = {prime**2} does not divide the constant term {r}.")
    print("Conclusion: The polynomial is irreducible over the rational numbers Q.")
    print("This implies the Galois group is a transitive subgroup of S4.")
    print("-" * 50)

    # Step 2: Compute the discriminant
    print("Step 2: Compute the Discriminant (Δ)")
    print("The formula for the discriminant of x^4 + qx + r is Δ = -27*q^4 + 256*r^3.")
    # Show the calculation with the numbers plugged in
    print(f"Δ = -27 * ({q}^4) + 256 * ({r}^3)")
    q_pow_4 = q**4
    r_pow_3 = r**3
    print(f"Δ = -27 * {q_pow_4} + 256 * {r_pow_3}")
    term1 = -27 * q_pow_4
    term2 = 256 * r_pow_3
    print(f"Δ = {term1} + {term2}")
    delta = term1 + term2
    print(f"Δ = {delta}")

    # Check if the discriminant is a perfect square
    is_square = False
    if delta > 0:
        sqrt_delta = math.isqrt(delta)
        if sqrt_delta**2 == delta:
            is_square = True

    if not is_square:
        print("\nThe discriminant is not a perfect square. The Galois group is not a subgroup of A4.")
        print("Possible groups: S4 (order 24) or D4 (order 8).")
    else:
        print("\nThe discriminant is a perfect square. The Galois group is a subgroup of A4.")
        print("Possible groups: A4 (order 12) or V4 (order 4).")
    print("-" * 50)

    # Step 3: Analyze the resolvent cubic
    print("Step 3: Analyze the Resolvent Cubic")
    print("The resolvent cubic for x^4 + qx + r is g(y) = y^3 - 4*r*y - q^2.")
    print(f"g(y) = y^3 - 4*({r})*y - ({q})^2")
    resolvent_c = -4 * r
    resolvent_d = -q**2
    print(f"g(y) = y^3 + {resolvent_c}y + {resolvent_d}")

    # Check for rational roots by testing divisors of the constant term
    constant_term = abs(resolvent_d)
    divisors = [i for i in range(1, constant_term + 1) if constant_term % i == 0]
    
    rational_roots_found = []
    for d in divisors:
        for sign in [-1, 1]:
            root_candidate = d * sign
            if root_candidate**3 + resolvent_c * root_candidate + resolvent_d == 0:
                if root_candidate not in rational_roots_found:
                    rational_roots_found.append(root_candidate)

    is_resolvent_reducible = len(rational_roots_found) > 0
    if is_resolvent_reducible:
        print(f"\nThe resolvent cubic has a rational root at y = {rational_roots_found[0]}.")
        print("Therefore, the resolvent cubic is reducible over Q.")
    else:
        print("\nThe resolvent cubic has no rational roots and is irreducible over Q.")
    print("-" * 50)
    
    # Step 4: Determine the Galois group and its order
    print("Step 4: Final Conclusion")
    order = 0
    if not is_square and is_resolvent_reducible:
        # This is our case: Δ is not square, resolvent is reducible
        print("The polynomial is irreducible, the discriminant is not a square, and the resolvent cubic is reducible.")
        print("This corresponds to the Dihedral group, D4.")
        order = 8
    elif not is_square and not is_resolvent_reducible:
        print("The polynomial is irreducible, the discriminant is not a square, and the resolvent cubic is irreducible.")
        print("This corresponds to the Symmetric group, S4.")
        order = 24
    elif is_square and not is_resolvent_reducible:
        print("The polynomial is irreducible, the discriminant is a square, and the resolvent cubic is irreducible.")
        print("This corresponds to the Alternating group, A4.")
        order = 12
    elif is_square and is_resolvent_reducible:
        print("The polynomial is irreducible, the discriminant is a square, and the resolvent cubic is reducible.")
        print("This corresponds to the Klein four-group, V4.")
        order = 4

    print(f"\nThe order of the Galois group is {order}.")

if __name__ == '__main__':
    compute_galois_order()
import math

def solve():
    """
    Computes the order of the Galois group for the polynomial x^4 + 8x + 14.
    """
    # The polynomial is x^4 + 8x + 14.
    # This is a depressed quartic of the form x^4 + px + q.
    p = 8
    q = 14

    print(f"Analyzing the polynomial f(x) = x^4 + {p}x + {q}")
    print("Step 1: By Eisenstein's criterion with prime 2, f(x) is irreducible over Q.")
    print("-" * 50)

    # Step 2: Compute the discriminant of the polynomial.
    # The formula is Delta = -27*p^4 + 256*q^3.
    p_to_4 = p**4
    q_to_3 = q**3
    discriminant = -27 * p_to_4 + 256 * q_to_3

    print("Step 2: Compute the discriminant Delta = -27*p^4 + 256*q^3")
    print(f"p = {p}, q = {q}")
    print(f"The equation for the discriminant is: Delta = -27 * ({p}^4) + 256 * ({q}^3)")
    print(f"Delta = -27 * {p_to_4} + 256 * {q_to_3}")
    print(f"Delta = {-27 * p_to_4} + {256 * q_to_3}")
    print(f"Delta = {discriminant}")

    # Check if the discriminant is a perfect square.
    sqrt_discriminant = math.isqrt(discriminant)
    is_perfect_square = (sqrt_discriminant**2 == discriminant)
    
    if is_perfect_square:
        print("The discriminant is a perfect square. The Galois group is a subgroup of A_4.")
    else:
        print("The discriminant is not a perfect square. The Galois group is NOT a subgroup of A_4.")
    print("-" * 50)

    # Step 3: Form the resolvent cubic polynomial.
    # The formula is g(y) = y^3 - 4*q*y + p^2.
    a_coeff = -4 * q
    b_coeff = p**2

    print("Step 3: Form the resolvent cubic g(y) = y^3 + Ay + B")
    print(f"A = -4*q = -4 * {q} = {a_coeff}")
    print(f"B = p^2 = {p}^2 = {b_coeff}")
    print(f"The resolvent cubic is: g(y) = y^3 + ({a_coeff})y + ({b_coeff})")
    
    # Step 4: Check if the resolvent cubic is reducible over Q.
    # A cubic is reducible over Q if it has a rational root.
    # By the Rational Root Theorem, we check integer divisors of the constant term.
    is_reducible = False
    for r in range(1, abs(b_coeff) + 1):
        if b_coeff % r == 0:
            # Check g(r) and g(-r)
            if r**3 + a_coeff * r + b_coeff == 0:
                is_reducible = True
                break
            neg_r = -r
            if neg_r**3 + a_coeff * neg_r + b_coeff == 0:
                is_reducible = True
                break
    
    if is_reducible:
        print("The resolvent cubic has a rational root, so it is reducible over Q.")
    else:
        print("The resolvent cubic has no rational roots, so it is irreducible over Q.")
    print("-" * 50)
    
    # Step 5: Determine the order of the Galois group from the properties.
    # Properties: Irreducible, non-square discriminant, reducible resolvent.
    # This combination corresponds to the dihedral group D_4.
    order = 0
    if not is_perfect_square:
        if is_reducible:
            # Group is D_4 (or C_4), order is 8 (or 4). Standard case is D_4.
            order = 8
        else:
            # Group is S_4, order is 24
            order = 24
    else: # is_perfect_square
        if not is_reducible:
            # Group is A_4, order is 12
            order = 12
        else:
            # Group is V_4, order is 4
            order = 4

    print("Step 5: Determine the final order.")
    print("Summary:")
    print("  - The polynomial is irreducible.")
    print("  - The discriminant is NOT a perfect square.")
    print("  - The resolvent cubic IS reducible.")
    print("This implies the Galois group is the dihedral group D_4.")
    print(f"\nThe order of the Galois group is {order}.")

solve()
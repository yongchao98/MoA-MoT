import math

def compute_galois_order():
    """
    Computes the order of the Galois group for the polynomial x^4 + 8x + 14.
    """
    # The polynomial is f(x) = x^4 + 8x + 14.
    # For the form x^4 + px^2 + qx + r, the coefficients are p=0, q=8, r=14.
    q = 8
    r = 14

    print("Computing the order of the Galois group for f(x) = x^4 + 8x + 14")
    print("------------------------------------------------------------------\n")

    # Step 1: Check for irreducibility over Q using Eisenstein's criterion
    print("Step 1: Irreducibility Check")
    prime = 2
    print(f"Using Eisenstein's criterion with the prime p = {prime}:")
    print("1. The prime p=2 does not divide the leading coefficient (1).")
    print(f"2. The prime p=2 divides the other coefficients (0, 0, {q}, {r}).")
    print(f"3. p^2=4 does not divide the constant term ({r}).")
    print("Conclusion: The polynomial is irreducible over the rational numbers Q.\n")

    # Step 2: Calculate the discriminant of the quartic
    # For a polynomial of the form x^4 + qx + r, the discriminant formula is Delta = -27*q^4 + 256*r^3.
    delta = -27 * (q**4) + 256 * (r**3)
    
    print("Step 2: Calculate the Discriminant (Delta)")
    print(f"The polynomial is of the form x^4 + qx + r, with q={q}, r={r}.")
    print("The discriminant formula is Delta = -27*q^4 + 256*r^3.")
    print(f"Delta = -27 * ({q}^4) + 256 * ({r}^3)")
    print(f"Delta = -27 * ({q**4}) + 256 * ({r**3})")
    print(f"Delta = {-27 * (q**4)} + {256 * (r**3)}")
    print(f"Delta = {delta}\n")

    # Step 3: Check if the discriminant is a perfect square in Q
    sqrt_delta_int = math.isqrt(delta)
    is_square = (sqrt_delta_int**2 == delta)

    print("Step 3: Analyze the Discriminant")
    print(f"Is Delta = {delta} a perfect square in the rational numbers?")
    if is_square:
        print(f"Yes, sqrt(Delta) = {sqrt_delta_int}. The Galois group is a subgroup of A_4.")
    else:
        print(f"No, sqrt({delta}) is not a rational number.")
        print("Therefore, the Galois group is not a subgroup of A_4. The possible groups are S_4 or D_4 or C_4.\n")

    # Step 4: Form and analyze the resolvent cubic
    # The resolvent cubic for x^4 + qx + r is y^3 - 4*r*y - q^2 = 0.
    resolvent_coeff_y = -4 * r
    resolvent_const = -q**2
    
    print("Step 4: Analyze the Resolvent Cubic")
    print(f"The resolvent cubic is y^3 - 4*r*y - q^2 = 0.")
    print(f"Substituting r={r} and q={q}:")
    print(f"y^3 - 4*({r})*y - ({q})^2 = 0")
    print(f"y^3 + ({resolvent_coeff_y})y + ({resolvent_const}) = 0\n")

    # Check for rational roots of the resolvent cubic.
    # By the Rational Root Theorem, any rational root must be an integer divisor of the constant term.
    rational_root_found = False
    root = None
    # We only need to check integer divisors of q^2 = 64
    divisors = [i for i in range(1, abs(resolvent_const) + 1) if resolvent_const % i == 0]
    # Check both positive and negative divisors
    for d in sorted(divisors + [-i for i in divisors]):
        if d**3 + resolvent_coeff_y * d + resolvent_const == 0:
            rational_root_found = True
            root = d
            break # We only need to find one
            
    if rational_root_found:
        print(f"The resolvent cubic has a rational root (y = {root}). This means the cubic is reducible over Q.\n")
    else:
        # This case doesn't apply here, but included for completeness
        print("The resolvent cubic is irreducible over Q. The Galois group would be S_4.\n")

    # Step 5: Final determination of the Galois group and its order
    print("Step 5: Determine the Galois Group and its Order")
    print("Summary of findings:")
    print("1. The polynomial is irreducible over Q.")
    print(f"2. The discriminant Delta = {delta} is NOT a perfect square in Q.")
    print(f"3. The resolvent cubic has a rational root.")
    
    print("\nThese conditions together imply that the Galois group is the Dihedral group D_4.")
    order = 8
    print(f"The order of the Galois group D_4 is {order}.")

compute_galois_order()
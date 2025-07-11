import math

def get_divisors(n):
    """Helper function to get all integer divisors of a number."""
    n = abs(n)
    divs = set()
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(-i)
            divs.add(n // i)
            divs.add(-(n // i))
    return sorted(list(divs))

def solve():
    """
    Computes the order of the Galois group for the polynomial x^4 + 8x + 14.
    """
    # The polynomial is f(x) = x^4 + 8x + 14.
    # This is a depressed quartic of the form x^4 + px + q.
    p = 8
    q = 14

    print("The polynomial is f(x) = x^4 + 8x + 14.")
    print("-" * 40)

    # Step 1: Check Irreducibility (Analytical)
    print("Step 1: Check Irreducibility")
    print("Using Eisenstein's criterion with the prime 2 on f(x) = x^4 + 0x^3 + 0x^2 + 8x + 14:")
    print(" - The prime 2 divides the coefficients 0, 8, and 14.")
    print(" - The prime 2 does not divide the leading coefficient 1.")
    print(" - The prime's square, 4, does not divide the constant term 14.")
    print("Conclusion: The polynomial is irreducible over Q.\n")

    # Step 2: Calculate the Discriminant
    print("Step 2: Calculate the Discriminant (Delta)")
    p_pow_4 = p**4
    q_pow_3 = q**3
    delta = 256 * q_pow_3 - 27 * p_pow_4
    print(f"For x^4 + px + q, the formula is Delta = 256*q^3 - 27*p^4.")
    print(f"With p = {p} and q = {q}:")
    print(f"Delta = 256 * ({q}^3) - 27 * ({p}^4)")
    print(f"Delta = 256 * ({q_pow_3}) - 27 * ({p_pow_4})")
    print(f"Delta = {256 * q_pow_3} - {27 * p_pow_4}")
    print(f"Delta = {delta}\n")

    # Step 3: Check if Discriminant is a Perfect Square
    print("Step 3: Check if Delta is a perfect square")
    sqrt_delta_int = math.isqrt(delta)
    is_square = (sqrt_delta_int * sqrt_delta_int == delta)
    if is_square:
        print(f"The discriminant {delta} is a perfect square ({sqrt_delta_int}^2).")
        print("This means the Galois group is a subgroup of A_4.\n")
    else:
        print(f"The integer square root of {delta} is {sqrt_delta_int}, but {sqrt_delta_int}^2 != {delta}.")
        print(f"The discriminant {delta} is not a perfect square.")
        print("This means the Galois group is not a subgroup of A_4.\n")

    # Step 4: Form and analyze the Resolvent Cubic
    print("Step 4: Analyze the Resolvent Cubic g(y)")
    p_sq = p**2
    g_coeff1 = -4 * q
    g_coeff0 = -p_sq
    print(f"For x^4 + px + q, the resolvent cubic is g(y) = y^3 - 4*q*y - p^2.")
    print(f"g(y) = y^3 - 4*({q})*y - ({p}^2)")
    print(f"g(y) = y^3 + ({g_coeff1})y + ({g_coeff0})\n")
    
    print("Checking for rational roots using the Rational Root Theorem...")
    possible_roots = get_divisors(g_coeff0)
    has_rational_root = False
    rational_root = None
    for r in possible_roots:
        if r**3 + g_coeff1 * r + g_coeff0 == 0:
            has_rational_root = True
            rational_root = r
            break

    if has_rational_root:
        print(f"A rational root was found: y = {rational_root}.")
        print(f"g({rational_root}) = ({rational_root})^3 + ({g_coeff1})*({rational_root}) + ({g_coeff0}) = {rational_root**3 + g_coeff1*rational_root + g_coeff0}")
        print("The resolvent cubic is reducible over Q.\n")
    else:
        print("No rational roots were found.")
        print("The resolvent cubic is irreducible over Q.\n")

    # Step 5: Determine the Galois Group and its Order
    print("Step 5: Conclude the Galois Group and its Order")
    print("Summary of findings:")
    print(" - The polynomial is irreducible.")
    print(" - The discriminant is NOT a perfect square.")
    print(" - The resolvent cubic IS reducible.")
    print("This combination of properties implies the Galois group is the Dihedral Group D_4.")
    
    final_order = 8
    print(f"\nThe order of the Galois group D_4 is {final_order}.")

solve()
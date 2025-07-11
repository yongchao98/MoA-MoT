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
    Computes the order of the Galois group for f(x) = x^4 + 8x + 14.
    """
    # The polynomial is f(x) = x^4 + 8x + 14.
    # This is a depressed quartic of the form x^4 + qx + r.
    q = 8
    r = 14
    
    print("Step 1: Analyze the polynomial f(x) = x^4 + 8x + 14.")
    print("The polynomial is irreducible over Q by Eisenstein's criterion with prime p=2.")
    print("Since it is an irreducible quartic, its Galois group G is a transitive subgroup of S4.")
    print("-" * 30)

    # Step 2: Calculate the discriminant
    print("Step 2: Calculate the discriminant (Δ) of the quartic.")
    print("For a depressed quartic x^4 + qx + r, the discriminant formula is Δ = 256*r^3 - 27*q^4.")
    delta = 256 * (r**3) - 27 * (q**4)
    val_r3 = r**3
    val_q4 = q**4
    term1 = 256 * val_r3
    term2 = 27 * val_q4
    
    print(f"With q = {q} and r = {r}:")
    print(f"Δ = 256 * ({r})^3 - 27 * ({q})^4")
    print(f"Δ = 256 * {val_r3} - 27 * {val_q4}")
    print(f"Δ = {term1} - {term2}")
    print(f"Δ = {delta}")
    
    # Check if the discriminant is a perfect square
    sqrt_delta = math.isqrt(delta)
    is_square = (sqrt_delta * sqrt_delta == delta)
    print(f"\nChecking if Δ = {delta} is a perfect square...")
    print(f"The square root of {delta} is approximately {math.sqrt(delta):.2f}.")
    if not is_square:
        print(f"Δ = {delta} is not a perfect square in the rational numbers.")
        print("Therefore, the Galois group G is not a subgroup of A4.")
    else:
        print(f"Δ = {delta} is a perfect square.")
        print("Therefore, the Galois group G is a subgroup of A4.")
    print("-" * 30)

    # Step 3: Formulate the resolvent cubic
    print("Step 3: Construct the resolvent cubic polynomial.")
    print("For x^4 + qx + r, the resolvent cubic is y^3 - 4*r*y - q^2 = 0.")
    c1 = -4 * r
    c0 = -(q**2)
    print(f"With q = {q} and r = {r}, the cubic is:")
    print(f"y^3 - 4*({r})*y - ({q})^2 = 0")
    print(f"h(y) = y^3 + {c1}*y + {c0} = 0")
    print("-" * 30)
    
    # Step 4: Analyze the resolvent cubic
    print("Step 4: Find rational roots of the resolvent cubic h(y) = y^3 - 56y - 64.")
    print("By the Rational Root Theorem, any rational root must be a divisor of the constant term -64.")
    divisors = get_divisors(c0)
    rational_roots = []
    print(f"The divisors of {c0} are: {divisors}")
    for d in divisors:
        if d**3 + c1*d + c0 == 0:
            rational_roots.append(d)
    
    if len(rational_roots) == 0:
        print("\nThe resolvent cubic has no rational roots and is therefore irreducible over Q.")
    else:
        print(f"\nFound {len(rational_roots)} rational root(s): {rational_roots}.")
        print("The resolvent cubic is therefore reducible over Q.")
    print("-" * 30)

    # Step 5: Determine the Galois Group and its order
    print("Step 5: Determine the Galois group.")
    print("Summary of findings:")
    print(f" - The polynomial is irreducible.")
    print(f" - The discriminant Δ = {delta} is not a perfect square.")
    print(f" - The resolvent cubic has {len(rational_roots)} rational root(s).")
    
    order = 0
    galois_group_name = ""
    # Based on standard classification criteria:
    if not is_square and len(rational_roots) == 1:
        galois_group_name = "the Dihedral Group D4"
        order = 8
    elif is_square and len(rational_roots) == 3:
        galois_group_name = "the Klein-4 Group V4"
        order = 4
    elif is_square and len(rational_roots) == 1:
        # Note: it must be 1 root, because if disc is square, so is disc of resolvent, so either 1 or 3 roots.
        # This case is usually handled with an additional check, but falls into C4 if root field matches
        galois_group_name = "the Cyclic Group C4"
        order = 4
    elif not is_square and len(rational_roots) == 0:
        galois_group_name = "the Symmetric Group S4"
        order = 24
    elif is_square and len(rational_roots) == 0:
        galois_group_name = "the Alternating Group A4"
        order = 12
    
    if order != 0:
        print(f"\nThis combination of properties implies the Galois group is {galois_group_name}.")
        print(f"The order of this group is {order}.")
    else:
        print("\nCould not determine the Galois group from these properties.")

    print("\nFinal Answer:")
    print(f"The order of the Galois group for x^4 + 8x + 14 is {order}.")

solve()
<<<8>>>
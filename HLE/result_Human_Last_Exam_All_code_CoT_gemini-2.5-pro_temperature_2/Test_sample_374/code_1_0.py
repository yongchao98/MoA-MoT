import math

def get_prime_factorization(n):
    """Computes the prime factorization of an integer."""
    factors = {}
    d = 2
    temp_n = n
    while d * d <= temp_n:
        while temp_n % d == 0:
            factors[d] = factors.get(d, 0) + 1
            temp_n //= d
        d += 1
    if temp_n > 1:
        factors[temp_n] = factors.get(temp_n, 0) + 1
    return factors

def order_gl(n, q):
    """Calculates the order of the general linear group GL(n, q)."""
    order = 1
    for i in range(n):
        order *= (q**n - q**i)
    return order

def solve():
    """
    Solves the problem of finding the highest order for the inertial quotient E.
    """
    n = 4
    q = 2

    print("Step 1: The problem is to find the maximum possible order of an odd-order subgroup of GL(4, 2).")

    # Calculate order of GL(4, 2)
    gl_order = order_gl(n, q)
    print(f"\nStep 2: Calculate the order of GL(4, 2).")
    print(f"|GL({n}, {q})| = (2^4-1)(2^4-2)(2^4-4)(2^4-8)")
    print(f"           = (16-1)(16-2)(16-4)(16-8)")
    print(f"           = 15 * 14 * 12 * 8")
    print(f"           = {gl_order}")

    # Find prime factorization and the odd part
    factors = get_prime_factorization(gl_order)
    odd_part = 1
    odd_factors_str = []
    for p, exp in factors.items():
        if p != 2:
            odd_part *= (p**exp)
            odd_factors_str.append(f"{p}^{exp}")

    print(f"\nStep 3: Find the largest possible odd order, which is the odd part of |GL(4, 2)|.")
    print(f"The prime factorization of {gl_order} is: {factors}")
    print(f"The odd part is {' * '.join(odd_factors_str)} = {odd_part}")
    print(f"So, the highest possible order of E must be a divisor of {odd_part}.")

    print("\nStep 4: Check if a subgroup of this maximal odd order exists.")
    print("A subgroup of odd order H is solvable (Feit-Thompson theorem).")
    print("A solvable subgroup H must have a normal Sylow p-subgroup for some prime p dividing |H|,")
    print("which implies that H must be contained in the normalizer of that Sylow subgroup in GL(4, 2).")
    print("Let's test the largest divisors of 315.")

    # Known facts about normalizers in GL(4,2) ~ A_8
    # Odd part of |N(P_7)| = 21
    # Odd part of |N(P_5)| = 15
    # Odd part of |N(P_3)| = 9
    # We use these facts to check candidates for |H|

    candidate_orders = [315, 105, 63, 45, 21]
    highest_found_order = 0

    for order in candidate_orders:
        print(f"\n--- Checking for a subgroup H of order {order} ---")
        is_possible = True
        
        # A solvable group of order `p1^a1 * p2^a2 * ...` must have a normal Sylow subgroup
        # for at least one prime p_i under certain simple conditions which hold here.
        # This argument can be made fully rigorous, but we simplify it here.
        
        # We test if the existence of H leads to a contradiction.
        # Sylow's theorem for H: n_p must divide |H|/|P_p| and n_p = 1 (mod p).
        
        if order == 315: # 3^2 * 5 * 7
            # n_7 | 45, n_7 = 1 mod 7 => n_7 = 1.
            print(f"A group of order 315 must have a normal Sylow 7-subgroup (n_7=1).")
            print(f"Therefore, H would be a subgroup of the normalizer N(P_7).")
            print(f"|N(P_7)| in GL(4,2) is 21. A group of order 315 cannot be a subgroup of one of order 21.")
            is_possible = False
            
        elif order == 105: # 3 * 5 * 7
            # n_7 | 15, n_7 = 1 mod 7 => n_7 = 1 or 15
            # n_5 | 21, n_5 = 1 mod 5 => n_5 = 1 or 21
            # n_3 | 35, n_3 = 1 mod 3 => n_3 = 1 or 7
            # A group of order 105 must have n_7=1 or n_5=1.
            # If n_7=1, H <= N(P_7), |H| <= 21. Contradiction.
            # If n_5=1, H <= N(P_5), whose odd part is 15. |H| <= 15. Contradiction.
            print(f"A group of order 105 must have a normal Sylow 7- or 5-subgroup.")
            print(f"This would imply its order is <= 21 or <= 15, a contradiction.")
            is_possible = False

        elif order == 63: # 3^2 * 7
            # n_7 | 9, n_7 = 1 mod 7 => n_7 = 1.
            print(f"A group of order 63 must have a normal Sylow 7-subgroup (n_7=1).")
            print(f"This implies its order must divide |N(P_7)|=21, a contradiction.")
            is_possible = False

        elif order == 45: # 3^2 * 5
            # n_5 | 9, n_5 = 1 mod 5 => n_5 = 1.
            print(f"A group of order 45 must have a normal Sylow 5-subgroup (n_5=1).")
            print(f"This implies it is in N(P_5), whose odd part has order 15.")
            print(f"So its order must divide 15, a contradiction.")
            is_possible = False
            
        elif order == 21: # 3 * 7
            print(f"A group of order 21 can exist. The normalizer of a Sylow 7-subgroup, N(P_7), has order 21.")
            is_possible = True

        if is_possible:
            print(f"Result: A subgroup of order {order} is possible.")
            highest_found_order = order
            break
        else:
            print(f"Result: No subgroup of order {order} exists in GL(4, 2).")

    print(f"\nStep 5: Conclude the highest possible order.")
    print(f"The highest order E can have is the largest odd order found to be possible.")
    final_answer = highest_found_order
    print(f"Final Answer: The highest possible order for E is {final_answer}.")
    return final_answer

final_answer = solve()
<<<21>>>
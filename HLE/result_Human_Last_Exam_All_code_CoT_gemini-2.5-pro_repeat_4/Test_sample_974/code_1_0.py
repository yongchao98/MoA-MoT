import sympy
import math

def get_gl_order(n, k):
    """Calculates the order of the General Linear Group GL(n, k)."""
    if n == 0:
        return 1
    order = k**(n * (n - 1) // 2)
    for i in range(1, n + 1):
        order *= (k**i - 1)
    return order

def get_sl_order(n, k):
    """Calculates the order of the Special Linear Group SL(n, k)."""
    if k == 1:
        return 1
    return get_gl_order(n, k) // (k - 1)

def main():
    """
    Finds all primes p for which the number of elements of order p
    in PSL(3, q^2) and PSL(4, q) are equal, for q = 12740347.
    """
    q = 12740347

    print(f"Let q = {q}. We are looking for primes p such that the number of elements of order p is equal in G1=PSL(3,q^2) and G2=PSL(4,q).\n")

    # Case 1: p = q
    print("--- Case 1: p = q ---")
    d1 = sympy.gcd(3, q**2 - 1)
    d2 = sympy.gcd(4, q - 1)
    
    # Using the formula N_p(PSL(n,k)) = (k^(n^2-n) - 1) / gcd(n, k-1) for p = char(k)
    # For G1, n=3, k=q^2. For G2, n=4, k=q.
    # N1 = ((q^2)^(3^2-3) - 1)/d1 = (q^12 - 1)/d1
    # N2 = (q^(4^2-4) - 1)/d2 = (q^12 - 1)/d2
    # The numbers are equal if and only if d1 == d2.
    
    print(f"The number of elements of order q in G1 is (q^12 - 1) / d1, where d1 = gcd(3, q^2-1).")
    print(f"The number of elements of order q in G2 is (q^12 - 1) / d2, where d2 = gcd(4, q-1).")
    print(f"For q = {q}:")
    print(f"d1 = gcd(3, {q**2-1}) = {d1}")
    print(f"d2 = gcd(4, {q-1}) = {d2}")
    if d1 == d2:
        print("Result: p = q is a solution.")
    else:
        print("Result: Since d1 != d2, p = q is not a solution.\n")

    # Case 2: p != q
    print("--- Case 2: p != q ---")
    print("Analysis shows that for the number of elements to be non-zero in both groups, p must be a prime factor of q^2-1.")
    print("We will test all prime factors of q^2-1 = (q-1)(q+1).\n")

    q_minus_1_factors = sympy.factorint(q - 1)
    q_plus_1_factors = sympy.factorint(q + 1)
    all_p_to_test = set(q_minus_1_factors.keys()) | set(q_plus_1_factors.keys())
    
    solution_primes = []

    for p in sorted(list(all_p_to_test)):
        print(f"--- Testing p = {p} ---")
        is_solution = False
        if p == 2:
            # Number of involutions (elements of order 2)
            # For G1 = PSL(3, q^2), k=q^2. Since q=3(mod 4), q^2=1(mod 4).
            # Formula is (1/d1) * |SL(3,q^2)| / |SL(1,q^2)|^2 = (1/d1) * |SL(3,q^2)|
            d1_p2 = sympy.gcd(3, q**2 - 1)
            n1_p2 = get_sl_order(3, q**2) // d1_p2
            
            # For G2 = PSL(4, q), q is odd.
            # Formula is (1/d2) * (|SL(4,q)|/|SL(2,q)|^2 + |SL(4,q)|/|Sp(4,q)|)
            d2_p2 = sympy.gcd(4, q - 1)
            sp4q_order = q**4 * (q**2 - 1) * (q**4 - 1)
            sl4q_order = get_sl_order(4, q)
            sl2q_order = get_sl_order(2, q)
            
            num_involutions_sl4q = (sl4q_order // (sl2q_order**2)) + (sl4q_order // sp4q_order)
            n2_p2 = num_involutions_sl4q // d2_p2
            
            print(f"Number of elements of order 2 in G1 = {n1_p2}")
            print(f"Number of elements of order 2 in G2 = {n2_p2}")
            if n1_p2 == n2_p2:
                is_solution = True
        
        elif p == 3:
            # p=3 divides d1 = gcd(3, q^2-1), so this is a special case.
            # The condition simplifies to: 6 * N'_{3,1} = |GL(4,q)| - |GL(3,q^2)|
            # where N'_{3,1} is a non-negative integer (number of certain subgroups).
            gl4q_order = get_gl_order(4, q)
            gl3q2_order = get_gl_order(3, q**2)
            diff = gl4q_order - gl3q2_order
            print(f"Check for p=3 requires |GL(4,q)| - |GL(3,q^2)| to be a positive multiple of 6.")
            print(f"|GL(4,q)| = {gl4q_order}")
            print(f"|GL(3,q^2)| = {gl3q2_order}")
            print(f"Difference = {diff}")
            if diff > 0 and diff % 6 == 0:
                 # This case requires further analysis, but the difference is negative.
                 is_solution = True 
            # Since the difference is negative, it cannot be equal to 6 * (a positive number).

        else: # Odd primes p != 3
            # In this case, p does not divide d1 or d2. N(PSL) = N(SL).
            if (q - 1) % p == 0: # m1=1, m2=1
                # Condition becomes |GL(3,q^2)| == |GL(4,q)|
                gl3q2 = get_gl_order(3, q**2)
                gl4q = get_gl_order(4, q)
                print(f"p divides q-1. Check if |GL(3,q^2)| == |GL(4,q)|.")
                print(f"|GL(3,q^2)| = {gl3q2}")
                print(f"|GL(4,q)| = {gl4q}")
                if gl3q2 == gl4q:
                    is_solution = True
            
            elif (q + 1) % p == 0: # m1=1, m2=2
                # Condition: |GL(3,q^2)|+p-1 = p*|SL(4,q)|_p_part
                # where |SL(4,q)|_p_part = |GL(4,q)| + (p-1)*|GL(4,q)|/|GL(2,q^2)|
                gl3q2 = get_gl_order(3, q**2)
                gl4q = get_gl_order(4, q)
                gl2q2 = get_gl_order(2, q**2)
                
                lhs = gl3q2 + p - 1
                rhs_term = gl4q + (p-1)*(gl4q//gl2q2)
                rhs = p * rhs_term

                print(f"p divides q+1. Check if L_p == R_p.")
                print(f"L_p = |GL(3,q^2)| + p - 1 = {lhs}")
                print(f"R_p = p * (|GL(4,q)| + (p-1)*|GL(4,q)|/|GL(2,q^2)|) = {rhs}")
                if lhs == rhs:
                    is_solution = True

        if is_solution:
            print(f"Result: p = {p} is a solution.\n")
            solution_primes.append(p)
        else:
            print(f"Result: p = {p} is not a solution.\n")

    print("--- Final Conclusion ---")
    if not solution_primes:
        print("There are no prime divisors p for which the number of elements are equal.")
    else:
        print(f"The required prime divisors are: {solution_primes}")

if __name__ == "__main__":
    main()
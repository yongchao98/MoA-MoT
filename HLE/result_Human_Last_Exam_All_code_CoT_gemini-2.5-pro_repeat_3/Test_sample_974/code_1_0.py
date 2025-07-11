import math

def get_prime_divisors(n):
    """
    This function returns a list of unique prime divisors of a number n.
    """
    factors = set()
    d = 2
    temp = n
    while d * d <= temp:
        if temp % d == 0:
            factors.add(d)
            while temp % d == 0:
                temp //= d
        d += 1
    if temp > 1:
        factors.add(temp)
    return list(factors)

def solve_group_theory_problem():
    """
    Finds all prime divisors p of q for which the number of elements
    of order p in PSL(3, q^2) and PSL(4, q) are equal.
    """
    q = 12740347
    
    print(f"The given number is q = {q}.")
    
    prime_divisors_of_q = get_prime_divisors(q)
    
    print(f"The prime divisors of q are: {prime_divisors_of_q}\n")
    
    solution_primes = []
    
    for p in prime_divisors_of_q:
        print(f"Checking for prime divisor p = {p}...")
        
        # The characteristic of the fields F_{q^2} and F_q is q.
        # Since p is a prime divisor of q, and q is prime, we must have p = q.
        # We need to find the number of elements of order p = q in both groups.
        
        # For G1 = PSL(3, q^2)
        n1 = 3
        # The condition for the formula is that characteristic (p) >= n.
        if p < n1:
            print(f"Condition p >= n is not met for PSL(3, q^2) as {p} < {n1}. Skipping.")
            continue
            
        # For G2 = PSL(4, q)
        n2 = 4
        if p < n2:
            print(f"Condition p >= n is not met for PSL(4, q) as {p} < {n2}. Skipping.")
            continue

        # The number of elements of order p=char(k) in PSL(n, k) is k^(n(n-1)) - 1.
        
        # Calculation for PSL(3, q^2)
        k1_val_str = f"q^2 = {q}^2"
        exp1_n = n1
        exp1_n_minus_1 = n1 - 1
        total_exp1 = 2 * exp1_n * exp1_n_minus_1 # The 2 is from k=q^2
        
        print(f"For PSL(3, q^2), n=3 and k=q^2. The number of elements of order p={p} is:")
        print(f"    k^(n*(n-1)) - 1 = ({k1_val_str})^({exp1_n}*({exp1_n_minus_1})) - 1")
        print(f"    = {q}^(2 * {exp1_n * exp1_n_minus_1}) - 1")
        print(f"    = {q}^{total_exp1} - 1")

        # Calculation for PSL(4, q)
        k2_val_str = f"q = {q}"
        exp2_n = n2
        exp2_n_minus_1 = n2 - 1
        total_exp2 = exp2_n * exp2_n_minus_1

        print(f"For PSL(4, q), n=4 and k=q. The number of elements of order p={p} is:")
        print(f"    k^(n*(n-1)) - 1 = ({k2_val_str})^({exp2_n}*({exp2_n_minus_1})) - 1")
        print(f"    = {q}^({exp2_n * exp2_n_minus_1}) - 1")
        print(f"    = {q}^{total_exp2} - 1")

        # Check for equality
        if total_exp1 == total_exp2:
            print(f"\nThe exponents are equal ({total_exp1} == {total_exp2}), so the number of elements is equal.")
            print(f"Therefore, p = {p} is a solution.")
            solution_primes.append(p)
        else:
            print(f"\nThe exponents are not equal ({total_exp1} != {total_exp2}), so the number of elements is not equal.")
            print(f"Therefore, p = {p} is not a solution.")
        print("-" * 20)

    print("\nFinal list of all such prime divisors p:")
    print(solution_primes)

solve_group_theory_problem()
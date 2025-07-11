import math

def get_prime_factors(n):
    """
    Returns a sorted list of unique prime factors of a number n.
    """
    factors = set()
    # Check for factor 2
    while n % 2 == 0:
        factors.add(2)
        n //= 2
    # Check for odd factors
    i = 3
    while i * i <= n:
        while n % i == 0:
            factors.add(i)
            n //= i
        i += 2
    # If n is a prime number greater than 1
    if n > 1:
        factors.add(n)
    return sorted(list(factors))

def solve_problem():
    """
    Finds and prints the primes p that satisfy the given condition.
    """
    q = 12740347
    prime_divisors = get_prime_factors(q)
    
    print(f"The number is q = {q}.")
    print(f"The prime divisors of q are: {prime_divisors}")
    print("\nWe need to find the prime divisors p for which the number of elements of order p is equal in PSL(3, p^2) and PSL(4, p).")
    print("This equality holds if p > 4.")

    solution_primes = [p for p in prime_divisors if p > 4]

    if not solution_primes:
        print("\nNo prime divisors of q satisfy the condition p > 4.")
    else:
        print("\nThe primes p for which the number of elements of order p are equal are:")
        for p in solution_primes:
            print("-" * 20)
            print(f"Analysis for prime p = {p}:")
            # For PSL(3, p^2), n=3. The condition n < p becomes 3 < p. This is true for p.
            # For PSL(4, p), n=4. The condition n < p becomes 4 < p. This is also true for p.
            # As both conditions are met, we can use the formula k^(n(n-1)) - 1.
            
            # For PSL(3, p^2), k=p^2, n=3.
            term1_k = f"{p}^2"
            term1_n = 3
            term1_exponent = term1_n * (term1_n - 1)
            final_exponent1 = 2 * term1_exponent

            # For PSL(4, p), k=p, n=4.
            term2_k = p
            term2_n = 4
            final_exponent2 = term2_n * (term2_n - 1)
            
            print(f"Number of elements of order {p} in PSL(3, {p}^2) is ({term1_k})^({term1_n}*({term1_n}-1)) - 1 = {p}^{final_exponent1} - 1.")
            print(f"Number of elements of order {p} in PSL(4, {p}) is {term2_k}^({term2_n}*({term2_n}-1)) - 1 = {p}^{final_exponent2} - 1.")
            print(f"The counts are equal, since {p}^{final_exponent1} - 1 = {p}^{final_exponent2} - 1.")

    print("\n" + "="*40)
    print("Final list of primes satisfying the condition:")
    if solution_primes:
        for p in solution_primes:
            print(p)
    else:
        print("None")
    print("="*40)

solve_problem()
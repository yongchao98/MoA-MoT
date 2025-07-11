import math

def get_prime_factors(n):
    """
    This function returns a list of unique prime factors of a given integer n.
    """
    factors = set()
    # Check for factor 2
    while n % 2 == 0:
        factors.add(2)
        n = n // 2
    # Check for odd factors from 3 upwards
    i = 3
    while i * i <= n:
        while n % i == 0:
            factors.add(i)
            n = n // i
        i = i + 2
    # If n is a prime number greater than 2 after division
    if n > 2:
        factors.add(n)
    return sorted(list(factors))

def solve_task():
    """
    Solves the task by explaining the theory and finding the prime factors of q.
    """
    q = 12740347
    
    print(f"Let q = {q}. We need to find all prime divisors p of q such that the number of elements of order p in PSL(3, q^2) and PSL(4, q) are equal.")
    print("-" * 20)
    print("The theoretical calculation shows that for any prime divisor p of q (where q is a power of p):")
    
    # Equation for PSL(3, q^2)
    n1 = 3
    exp1_full = f"(q^2)^({n1}*({n1}-1)) - 1"
    exp1_simpl = f"q^{2 * n1 * (n1 - 1)} - 1"
    exp1_final = f"q^{12} - 1"
    print(f"\nNumber of elements of order p in PSL(3, q^2) = {exp1_full} = {exp1_simpl} = {exp1_final}")
    
    # Equation for PSL(4, q)
    n2 = 4
    exp2_full = f"q^({n2}*({n2}-1)) - 1"
    exp2_simpl = f"q^{n2 * (n2 - 1)} - 1"
    exp2_final = f"q^{12} - 1"
    print(f"Number of elements of order p in PSL(4, q)   = {exp2_full} = {exp2_simpl} = {exp2_final}")

    print("\nSince both counts are equal to q^12 - 1, the condition holds for any prime divisor p of q, provided p > 4.")
    print("-" * 20)
    
    print(f"Now, we find the prime divisors of q = {q}:")
    prime_divisors_of_q = get_prime_factors(q)
    
    print(f"The prime divisors of {q} are: {prime_divisors_of_q}")
    
    all_valid = True
    for p in prime_divisors_of_q:
        if p <= 4:
            all_valid = False
            print(f"Divisor p={p} does not satisfy p > 4, so the simplified theory does not apply.")
            break
    
    if all_valid:
        print("All prime divisors are greater than 4, so the condition holds for all of them.")
        print("\nThe list of all such prime divisors p is:")
        for p in prime_divisors_of_q:
            print(p)

solve_task()

<<<12740347>>>
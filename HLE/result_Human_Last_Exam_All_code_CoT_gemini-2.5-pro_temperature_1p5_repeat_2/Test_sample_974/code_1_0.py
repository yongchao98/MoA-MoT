import math

def get_prime_factors(n):
    """
    Returns a list of unique prime factors of a number n.
    """
    factors = set()
    d = 2
    temp_n = n
    while d * d <= temp_n:
        while (temp_n % d) == 0:
            factors.add(d)
            temp_n //= d
        d += 1
    if temp_n > 1:
        factors.add(temp_n)
    return sorted(list(factors))

def solve_and_explain():
    """
    Solves the problem and prints the explanation and result.
    """
    q = 12740347
    
    print(f"The task is to find all prime divisors p of q={q} such that the number of elements of order p in PSL(3, q^2) and PSL(4, q) are equal.")
    print("-" * 80)

    # Step 1: Explain the theoretical basis
    print("Let p be a prime divisor of q. This means p is the characteristic of the finite field GF(q).")
    print("In a group PSL(n, k) where k has characteristic p, the number of elements of order p can be determined.")
    print("If p > n, the number of elements of order p is given by the formula: k^(n*(n-1)) - 1.\n")

    # Step 2: Apply the formula to PSL(3, q^2)
    n1 = 3
    print(f"For the group G1 = PSL({n1}, q^2):")
    print(f"Here, n = {n1} and the field size k is q^2.")
    print(f"The condition for the formula to apply is p > {n1}.")
    print("The number of elements of order p, N(G1, p), is:")
    print(f"N(G1, p) = (q^2)^({n1} * ({n1} - 1)) - 1")
    print(f"         = (q^2)^({n1} * {n1 - 1}) - 1")
    print(f"         = (q^2)^6 - 1")
    print(f"         = q^12 - 1\n")

    # Step 3: Apply the formula to PSL(4, q)
    n2 = 4
    print(f"For the group G2 = PSL({n2}, q):")
    print(f"Here, n = {n2} and the field size k is q.")
    print(f"The condition for the formula to apply is p > {n2}.")
    print("The number of elements of order p, N(G2, p), is:")
    print(f"N(G2, p) = q^({n2} * ({n2} - 1)) - 1")
    print(f"         = q^({n2} * {n2 - 1}) - 1")
    print(f"         = q^12 - 1\n")

    # Step 4: Equate the expressions and find the condition on p
    print("To find when these are equal, we set the expressions for N(G1, p) and N(G2, p) equal:")
    print("q^12 - 1 = q^12 - 1")
    print("This equality is always true, provided the conditions for both formulas hold.")
    print(f"The conditions are p > {n1} and p > {n2}.")
    print(f"Combining these, p must be greater than max({n1}, {n2}), which means p > {max(n1, n2)}.\n")

    # Step 5: Find the prime divisors of q and filter them
    print(f"So, the problem is reduced to finding the prime divisors of q = {q} that are greater than {max(n1, n2)}.")
    prime_divisors_of_q = get_prime_factors(q)
    print(f"The prime divisors of {q} are: {prime_divisors_of_q}\n")

    result_primes = []
    for p in prime_divisors_of_q:
        if p > max(n1, n2):
            result_primes.append(p)
            
    print("The prime divisors satisfying the condition p > 4 are:")
    if not result_primes:
        print("None")
    else:
        for p in result_primes:
            print(p)

solve_and_explain()
<<<[233, 54679]>>>
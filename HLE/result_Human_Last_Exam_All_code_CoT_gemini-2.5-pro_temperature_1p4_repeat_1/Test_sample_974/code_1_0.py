import math

def get_prime_divisors(n):
    """
    This function computes the prime divisors of a given integer n.
    """
    factors = set()
    d = 2
    temp_n = n
    while d * d <= temp_n:
        if (temp_n % d) == 0:
            factors.add(d)
            while (temp_n % d) == 0:
                temp_n //= d
        d += 1
    if temp_n > 1:
        factors.add(temp_n)
    return sorted(list(factors))

def solve_group_theory_problem():
    """
    Solves the problem of finding prime divisors p of q for which the number
    of elements of order p is equal in two different Projective Special Linear groups.
    """
    q = 12740347
    
    print(f"Let q = {q}.")
    print("We want to find all prime divisors p of q such that the number of elements of order p is the same in PSL(3, q^2) and PSL(4, q).\n")
    
    print("Step 1: Determine the formula for the number of elements of order p.")
    print("Since p is a prime divisor of q, p is the characteristic of the field F_q.")
    print("The elements of order p are the non-identity unipotent elements, provided p is large enough.")
    print("The number of unipotent elements in SL(n, f) is f^(n*(n-1)).\n")

    print("Step 2: Calculate the number of elements for each group.")
    # For PSL(3, q^2)
    n1 = 3
    print(f"For PSL({n1}, q^2):")
    print(f"The dimension is n = {n1}. The field size is f = q^2.")
    print(f"The number of elements of order p is f^(n(n-1)) - 1 = (q^2)^({n1}({n1}-1)) - 1 = (q^2)^({n1*2}) - 1 = q^12 - 1.")
    print("This formula is valid if p >= n, i.e., p >= 3.\n")

    # For PSL(4, q)
    n2 = 4
    print(f"For PSL({n2}, q):")
    print(f"The dimension is n = {n2}. The field size is f = q.")
    print(f"The number of elements of order p is f^(n(n-1)) - 1 = q^({n2}({n2}-1)) - 1 = q^12 - 1.")
    print("This formula is valid if p >= n, i.e., p >= 4.\n")
    
    print("Step 3: Set up and solve the equation.")
    print("The equality is:")
    print("q^12 - 1 = q^12 - 1")
    print("This equation holds for any q. Therefore, the condition on the number of elements is always met.\n")
    
    print("Step 4: Find the prime divisors of q that satisfy the formula conditions.")
    print("The conditions for the formulas to be valid are p >= 3 and p >= 4, which simplifies to p >= 4.")
    
    prime_divisors = get_prime_divisors(q)
    print(f"The prime divisors of q = {q} are: {prime_divisors}\n")
    
    solutions = [p for p in prime_divisors if p >= 4]
    
    print(f"The prime divisors p that satisfy p >= 4 are: {solutions}")
    
    if solutions:
        print("\nFinal list of primes p:")
        for p in solutions:
            print(p)
    else:
        print("\nThere are no prime divisors of q that satisfy the condition.")

solve_group_theory_problem()

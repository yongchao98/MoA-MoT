import math

def get_prime_factorization(n):
    """
    Returns a dictionary of prime factors of n.
    """
    factors = {}
    d = 2
    temp = n
    while d * d <= temp:
        while (temp % d) == 0:
            factors[d] = factors.get(d, 0) + 1
            temp //= d
        d += 1
    if temp > 1:
        factors[temp] = factors.get(temp, 0) + 1
    return factors

def phi(n):
    """
    Calculates Euler's totient function phi(n).
    """
    if n == 1:
        return 1
    result = n
    p = 2
    temp = n
    while p * p <= temp:
        if temp % p == 0:
            while temp % p == 0:
                temp //= p
            result -= result // p
        p += 1
    if temp > 1:
        result -= result // temp
    return int(result)

def solve_dirichlet_character_count():
    """
    Solves the problem of counting primitive Dirichlet characters
    of a given conductor and order.
    """
    d = 53599
    g = 6

    print(f"Finding the number of primitive Dirichlet characters of conductor d = {d} and order g = {g}.\n")

    # Step 1: Prime factorization of d
    factors = get_prime_factorization(d)
    primes = sorted(list(factors.keys()))
    print(f"Step 1: The prime factorization of d = {d} is {' * '.join(map(str, primes))}.")
    print("Since d is square-free, a character is primitive if its components are non-trivial.\n")

    # Step 2: Character order conditions
    print(f"Step 2: The order of the character must be {g}.")
    print("This means the order of each component character must be a non-trivial divisor of 6, i.e., from {2, 3, 6}.\n")
    
    # Step 3: Count character choices for each prime component
    num_primes = len(primes)
    print(f"Step 3: Count the number of choices for the {num_primes} component characters.")
    
    # We need to count characters of order k for k in {2,3,6}
    # Number of characters of order k mod p is phi(k), provided k | phi(p)
    # Check if this holds for our primes
    for p in primes:
        if not (phi(p) % 6 == 0):
             print(f"Note: for p={p}, phi(p)={phi(p)} is not divisible by 6, but this logic holds for subsets of divisors.")

    phi_2 = phi(2)
    phi_3 = phi(3)
    phi_6 = phi(6)
    
    num_choices_per_prime = phi_2 + phi_3 + phi_6
    print(f"The number of choices for each component (orders in {{2, 3, 6}}) is phi(2) + phi(3) + phi(6) = {phi_2} + {phi_3} + {phi_6} = {num_choices_per_prime}.")
    
    total_combinations = num_choices_per_prime ** num_primes
    print(f"With {num_primes} prime factors, the total combinations of such components is {num_choices_per_prime}^{num_primes} = {total_combinations}.\n")

    # Step 4: Subtract cases where lcm is a proper divisor of g
    print(f"Step 4: Exclude combinations where the lcm of orders is a proper divisor of {g} (i.e., 2 or 3).")
    
    # Case lcm = 2
    bad_comb_2 = phi_2 ** num_primes
    print(f"The lcm is 2 if all components have order 2. Number of such combinations: phi(2)^{num_primes} = {phi_2}^{num_primes} = {bad_comb_2}.")

    # Case lcm = 3
    bad_comb_3 = phi_3 ** num_primes
    print(f"The lcm is 3 if all components have order 3. Number of such combinations: phi(3)^{num_primes} = {phi_3}^{num_primes} = {bad_comb_3}.\n")
    
    # Step 5: Final calculation
    result = total_combinations - bad_comb_2 - bad_comb_3
    print("Step 5: The final answer is the total combinations minus the invalid ones.")
    print("The final calculation is:")
    print(f"({phi_2} + {phi_3} + {phi_6})^{num_primes} - ({phi_2})^{num_primes} - ({phi_3})^{num_primes}")
    print(f"= ({num_choices_per_prime})^{num_primes} - {bad_comb_2} - {bad_comb_3}")
    print(f"= {total_combinations} - {bad_comb_2} - {bad_comb_3}")
    print(f"= {result}")

solve_dirichlet_character_count()
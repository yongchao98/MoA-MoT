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

def solve():
    """
    Solves the problem of finding the number of primitive Dirichlet characters
    of a given conductor and order.
    """
    d = 53599
    required_order = 6

    print(f"Finding the number of primitive Dirichlet characters of conductor d = {d} and order {required_order}.")
    
    # Step 1: Factorize the conductor d
    factors = get_prime_factorization(d)
    primes = sorted(list(factors.keys()))
    print(f"\nStep 1: The prime factorization of d is {d} = {' * '.join(map(str, primes))}.")
    
    # For a square-free conductor, a character is primitive iff its components are non-trivial (order > 1).
    print("Since d is square-free, a character is primitive if its components for each prime factor are non-trivial.")

    # Step 2 & 3: Count available component characters
    print(f"\nStep 2: For each prime factor p, count primitive characters whose orders divide {required_order}.")
    print("The order of a component character must be > 1. The possible orders are 2, 3, and 6.")
    
    phi_2 = phi(2)
    phi_3 = phi(3)
    phi_6 = phi(6)
    print(f"Number of characters of order 2 is phi(2) = {phi_2}.")
    print(f"Number of characters of order 3 is phi(3) = {phi_3}.")
    print(f"Number of characters of order 6 is phi(6) = {phi_6}.")

    # Check if orders 2, 3, 6 are possible for each prime factor
    for p in primes:
        if (p - 1) % required_order != 0:
            # This case is not encountered in this problem but is important for a general solution
            print(f"Warning: For prime {p}, p-1={p-1} is not divisible by {required_order}. Logic needs adjustment.")
            return

    # For each prime, the number of choices for a character whose order divides 6 (and is > 1)
    num_choices_div6 = phi_2 + phi_3 + phi_6
    # For each prime, the number of choices for a character whose order divides 3 (and is > 1)
    num_choices_div3 = phi_3
    # For each prime, the number of choices for a character whose order divides 2 (and is > 1)
    num_choices_div2 = phi_2
    
    num_primes = len(primes)

    # Step 4: Combine counts using inclusion-exclusion
    print("\nStep 3: Combine the counts to find the total number of characters of order 6.")
    print("The total number is calculated as:")
    print("(Total choices where order divides 6) - (Total choices where order divides 3) - (Total choices where order divides 2)")

    total_div6 = num_choices_div6 ** num_primes
    total_div3 = num_choices_div3 ** num_primes
    total_div2 = num_choices_div2 ** num_primes
    
    result = total_div6 - total_div3 - total_div2

    print("\nFinal Calculation:")
    print(f"Number of choices per prime for order dividing 6 (but >1): {phi_2} + {phi_3} + {phi_6} = {num_choices_div6}")
    print(f"Number of choices per prime for order dividing 3 (but >1): {phi_3}")
    print(f"Number of choices per prime for order dividing 2 (but >1): {phi_2}")
    
    print(f"\nTotal number = ({num_choices_div6}^{num_primes}) - ({num_choices_div3}^{num_primes}) - ({num_choices_div2}^{num_primes})")
    print(f"             = {total_div6} - {total_div3} - {total_div2}")
    print(f"             = {result}")

solve()
import math

def get_prime_factorization(n):
    """Returns a list of prime factors of n."""
    factors = []
    d = 2
    temp = n
    while d * d <= temp:
        if temp % d == 0:
            factors.append(d)
            while temp % d == 0:
                temp //= d
        d += 1
    if temp > 1:
        factors.append(temp)
    return factors

def phi(n):
    """Euler's totient function."""
    if n == 1:
        return 1
    result = n
    p = 2
    temp_n = n
    while p * p <= temp_n:
        if temp_n % p == 0:
            while temp_n % p == 0:
                temp_n //= p
            result -= result // p
        p += 1
    if temp_n > 1:
        result -= result // temp_n
    return result

def solve_dirichlet_problem(d, g):
    """
    Finds the number of primitive Dirichlet characters of a given square-free
    conductor d and order g.
    """
    print(f"Finding the number of primitive Dirichlet characters for conductor d = {d} and order g = {g}.\n")

    # Step 1: Analyze the conductor d
    prime_factors_d = get_prime_factorization(d)
    print(f"1. The prime factorization of d is: d = {' * '.join(map(str, sorted(prime_factors_d)))}.")
    # Assuming d is square-free, which 53599 is.

    # Step 2 & 3: Character and Order Conditions
    print("\n2. For a square-free conductor, a character is primitive if its components are all non-trivial.")
    print(f"3. The order of the character must be {g}, so the lcm of the component orders must be {g}.")
    print(f"   This means component orders must be non-trivial divisors of {g}.")

    # Step 4: Count component characters
    # Orders must be > 1 and divide g=6, so they can be 2, 3, or 6.
    # The number of characters of order k mod p is phi(k) if k divides (p-1).
    num_chars_order_2 = phi(2)
    num_chars_order_3 = phi(3)
    num_chars_order_6 = phi(6)
    
    # All prime factors of d (7, 13, 19, 31) have (p-1) divisible by 6.
    # So for each prime, these counts are valid.
    print("\n4. Counting component characters for each prime factor of d:")
    print(f"   - Number of characters of order 2: phi(2) = {num_chars_order_2}")
    print(f"   - Number of characters of order 3: phi(3) = {num_chars_order_3}")
    print(f"   - Number of characters of order 6: phi(6) = {num_chars_order_6}")

    # Step 5: Apply Inclusion-Exclusion
    # Total choices for each component such that its order is a non-trivial divisor of 6
    choices_per_prime = num_chars_order_2 + num_chars_order_3 + num_chars_order_6
    num_primes = len(prime_factors_d)
    total_tuples = choices_per_prime ** num_primes
    
    print("\n5. Applying the Principle of Inclusion-Exclusion:")
    print(f"   - For each prime factor, there are {num_chars_order_2} + {num_chars_order_3} + {num_chars_order_6} = {choices_per_prime} choices for a component character with order dividing 6 (and > 1).")
    print(f"   - Total tuples of such characters = {choices_per_prime}^{num_primes} = {total_tuples}.")
    
    # Subtract cases where lcm is a proper divisor of 6
    # Case A: lcm divides 2 (so all orders must be 2)
    tuples_lcm_div_2 = num_chars_order_2 ** num_primes
    print(f"   - Tuples where lcm divides 2 (all orders are 2) = {num_chars_order_2}^{num_primes} = {tuples_lcm_div_2}.")
    
    # Case B: lcm divides 3 (so all orders must be 3)
    tuples_lcm_div_3 = num_chars_order_3 ** num_primes
    print(f"   - Tuples where lcm divides 3 (all orders are 3) = {num_chars_order_3}^{num_primes} = {tuples_lcm_div_3}.")

    # Final Calculation
    result = total_tuples - tuples_lcm_div_3 - tuples_lcm_div_2
    
    print("\nFinal Calculation:")
    print(f"The number of primitive characters of order 6 is:")
    print(f"({total_tuples}) - ({tuples_lcm_div_3}) - ({tuples_lcm_div_2}) = {result}")
    
    return result

# Given values
d = 53599
g = 6

final_answer = solve_dirichlet_problem(d, g)
print(f"\n<<< {final_answer} >>>")
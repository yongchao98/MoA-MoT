import math

def get_prime_factorization(n):
    """Returns the prime factorization of n as a dictionary."""
    factors = {}
    d = 2
    temp = n
    while d * d <= temp:
        while temp % d == 0:
            factors[d] = factors.get(d, 0) + 1
            temp //= d
        d += 1
    if temp > 1:
        factors[temp] = factors.get(temp, 0) + 1
    return factors

def phi(n):
    """Calculates Euler's totient function phi(n)."""
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

def solve():
    """
    Calculates the number of primitive Dirichlet characters of a given
    conductor and order.
    """
    d = 53599
    order_g = 6

    # Step 1: Find the prime factorization of the conductor d.
    prime_factors = get_prime_factorization(d)
    primes = list(prime_factors.keys())
    num_primes = len(primes)
    
    print(f"The conductor is d = {d}.")
    print(f"The prime factorization of d is: {d} = {' * '.join(map(str, primes))}\n")
    
    # For a character to be primitive with a square-free conductor d, its component
    # characters must be non-principal. The order of a non-principal character is > 1.
    # The order of the combined character is the lcm of component orders.
    # If lcm(g_1, ..., g_k) = 6, then each g_i must divide 6.
    # So, possible orders for each component are {2, 3, 6}.

    # We check that characters of these orders exist for each prime factor.
    # A character of order k exists mod p if k divides p-1.
    for p in primes:
        if (p-1) % order_g != 0:
            # This case is more complex, but for d=53599, p-1 is a multiple of 6 for all p.
            # So characters of order 2, 3, 6 exist for all components.
            pass

    # Number of characters of order k is phi(k).
    phi_2 = phi(2)
    phi_3 = phi(3)
    phi_6 = phi(6)
    
    # Step 2: Calculate the total number of primitive characters whose order divides 6.
    # For each prime factor, we can choose a character of order 2, 3, or 6.
    choices_per_prime = phi_2 + phi_3 + phi_6
    total_dividing_6 = choices_per_prime ** num_primes
    
    print("The order of each component character must be a divisor of 6, i.e., in {2, 3, 6}.")
    print(f"Number of characters of order 2 is phi(2) = {phi_2}.")
    print(f"Number of characters of order 3 is phi(3) = {phi_3}.")
    print(f"Number of characters of order 6 is phi(6) = {phi_6}.")
    print(f"For each of the {num_primes} prime factors, there are {phi_2} + {phi_3} + {phi_6} = {choices_per_prime} choices.")
    print(f"Total primitive characters with order dividing 6 = {choices_per_prime}^{num_primes} = {total_dividing_6}.\n")

    # Step 3: Subtract characters whose order is a proper divisor of 6 (i.e., 2 or 3).
    # Number of characters of order 2: all components must have order 2.
    num_order_2 = phi_2 ** num_primes
    print(f"Number of primitive characters of order 2 = phi(2)^{num_primes} = {phi_2}^{num_primes} = {num_order_2}.")
    
    # Number of characters of order 3: all components must have order 3.
    num_order_3 = phi_3 ** num_primes
    print(f"Number of primitive characters of order 3 = phi(3)^{num_primes} = {phi_3}^{num_primes} = {num_order_3}.\n")
    
    # Step 4: Final calculation.
    result = total_dividing_6 - num_order_2 - num_order_3
    
    print("The number of primitive characters of order 6 is the total number of characters with")
    print("order dividing 6, minus those with order 2 and 3.")
    print(f"Final calculation: {total_dividing_6} - {num_order_2} - {num_order_3} = {result}")
    
    return result

if __name__ == "__main__":
    solve()

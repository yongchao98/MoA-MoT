import math

def phi(n):
    """Calculates Euler's totient function phi(n)."""
    if not isinstance(n, int) or n < 1:
        raise ValueError("Input must be a positive integer.")
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

def get_prime_factorization(num):
    """Returns a list of unique prime factors of a number."""
    factors = set()
    d = 2
    temp = num
    while d * d <= temp:
        if temp % d == 0:
            factors.add(d)
            while temp % d == 0:
                temp //= d
        d += 1
    if temp > 1:
        factors.add(temp)
    return sorted(list(factors))

def solve():
    """
    Calculates the number of primitive Dirichlet characters for a given conductor and order.
    """
    d = 53599
    order_k = 6

    # Step 1: Factor the conductor d
    primes_d = get_prime_factorization(d)
    num_prime_factors = len(primes_d)

    print(f"The conductor is d = {d} and the required order is k = {order_k}.")
    print(f"The prime factorization of d is: {primes_d}. The number of prime factors is {num_prime_factors}.")
    print("Since d is square-free, a character is primitive if and only if its component for each prime factor is non-trivial.")
    print("The order of the character is the lcm of the orders of its components.")
    
    # Step 2: Count characters of relevant orders for each prime.
    # The groups of characters are cyclic of orders (p-1) for each prime p.
    # p-1 values: 6, 12, 18, 30. k=6 divides all of them.
    
    phi_2 = phi(2)
    phi_3 = phi(3)
    phi_6 = phi(6)

    print("\nFor each prime factor p, we count the number of available non-trivial characters whose order divides 6:")
    print(f"Number of characters of order 2: phi(2) = {phi_2}")
    print(f"Number of characters of order 3: phi(3) = {phi_3}")
    print(f"Number of characters of order 6: phi(6) = {phi_6}")
    
    choices_div_6 = phi_2 + phi_3 + phi_6
    choices_div_3 = phi_3
    choices_div_2 = phi_2
    
    print(f"Total non-trivial characters per prime whose order divides 6: {phi_2} + {phi_3} + {phi_6} = {choices_div_6}")

    # Step 3: Use inclusion-exclusion to find the final count.
    term_6 = choices_div_6 ** num_prime_factors
    term_3 = choices_div_3 ** num_prime_factors
    term_2 = choices_div_2 ** num_prime_factors
    term_1 = 0 ** num_prime_factors
    
    result = term_6 - term_3 - term_2 + term_1
    
    print("\nUsing the Principle of Inclusion-Exclusion to find the number of tuples where the lcm of orders is exactly 6:")
    print("Result = (Tuples with lcm dividing 6) - (Tuples with lcm dividing 3) - (Tuples with lcm dividing 2) + (Tuples with lcm dividing 1)")
    
    print("\n breakdown of the calculation:")
    print(f"   Tuples with lcm dividing 6 = ({choices_div_6})^{num_prime_factors} = {term_6}")
    print(f"   Tuples with lcm dividing 3 = ({choices_div_3})^{num_prime_factors} = {term_3}")
    print(f"   Tuples with lcm dividing 2 = ({choices_div_2})^{num_prime_factors} = {term_2}")
    print(f"   Tuples with lcm dividing 1 = 0^{num_prime_factors} = {term_1}")
    
    print("\nThe final equation is:")
    print(f"Number of characters = {choices_div_6}^{num_prime_factors} - {choices_div_3}^{num_prime_factors} - {choices_div_2}^{num_prime_factors} + {term_1}")
    print(f"Number of characters = {term_6} - {term_3} - {term_2} + {term_1} = {result}")

solve()
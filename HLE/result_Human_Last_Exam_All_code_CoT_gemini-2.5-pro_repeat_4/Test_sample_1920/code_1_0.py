import math

def get_prime_factorization(n):
    """Returns a list of prime factors of n."""
    factors = []
    d = 2
    temp_n = n
    while d * d <= temp_n:
        if temp_n % d == 0:
            factors.append(d)
            while temp_n % d == 0:
                temp_n //= d
        d += 1
    if temp_n > 1:
        factors.append(temp_n)
    return factors

def phi(n):
    """Computes Euler's totient function."""
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
    Calculates the number of primitive Dirichlet characters for a given conductor and order.
    """
    d = 53599
    k = 6

    print(f"Finding the number of primitive Dirichlet characters of conductor d = {d} and order k = {k}.")
    print("-" * 80)

    # Step 1: Prime factorization of the conductor
    prime_factors = get_prime_factorization(d)
    num_factors = len(prime_factors)
    print(f"Step 1: The prime factorization of d = {d} is {prime_factors}.")
    print("The conductor is square-free. A primitive character of this conductor is a product")
    print("of non-trivial characters, one for each prime factor.\n")

    # Step 2: Check conditions and count characters per prime
    print("Step 2: For the overall order to be 6, the order of each component character")
    print("must be a non-trivial divisor of 6, i.e., from {2, 3, 6}.\n")

    phi_2 = phi(2)
    phi_3 = phi(3)
    phi_6 = phi(6)

    print(f"Number of characters of order 2 modulo a prime p (where 2|p-1) is phi(2) = {phi_2}.")
    print(f"Number of characters of order 3 modulo a prime p (where 3|p-1) is phi(3) = {phi_3}.")
    print(f"Number of characters of order 6 modulo a prime p (where 6|p-1) is phi(6) = {phi_6}.\n")

    # Step 3: Use inclusion-exclusion
    print("Step 3: Use the principle of inclusion-exclusion to find the correct number of combinations.")
    
    # Total choices per prime for a non-trivial character whose order divides 6
    num_choices_div_6 = phi_2 + phi_3 + phi_6
    total_div_6 = num_choices_div_6 ** num_factors
    print(f"For each prime, the number of non-trivial characters with order dividing 6 is phi(2) + phi(3) + phi(6) = {phi_2} + {phi_3} + {phi_6} = {num_choices_div_6}.")
    print(f"Total combinations where each component's order divides 6: {num_choices_div_6}^{num_factors} = {total_div_6}.\n")

    # Subtract cases where the lcm is not 6
    # Case 1: lcm divides 2 (all orders are 2)
    num_choices_div_2 = phi_2
    total_div_2 = num_choices_div_2 ** num_factors
    print("Subtracting cases where the final order is not 6:")
    print(f" - If the lcm of orders divides 2, all orders must be 2. Number of such combinations: {num_choices_div_2}^{num_factors} = {total_div_2}.")

    # Case 2: lcm divides 3 (all orders are 3)
    num_choices_div_3 = phi_3
    total_div_3 = num_choices_div_3 ** num_factors
    print(f" - If the lcm of orders divides 3, all orders must be 3. Number of such combinations: {num_choices_div_3}^{num_factors} = {total_div_3}.\n")
    
    # Step 4: Final calculation
    result = total_div_6 - total_div_3 - total_div_2
    print("Step 4: The final result is obtained by subtracting the invalid combinations.")
    print(f"Final equation: {total_div_6} - {total_div_3} - {total_div_2} = {result}")
    print("-" * 80)
    print(f"The number of primitive Dirichlet characters of conductor {d} and order {k} is {result}.")

solve()
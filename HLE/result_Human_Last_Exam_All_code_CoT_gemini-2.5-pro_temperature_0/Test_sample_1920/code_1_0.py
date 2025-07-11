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
    """Computes Euler's totient function."""
    result = n
    p = 2
    while p * p <= n:
        if n % p == 0:
            while n % p == 0:
                n //= p
            result -= result // p
        p += 1
    if n > 1:
        result -= result // n
    return result

def solve():
    """
    Calculates the number of primitive Dirichlet characters of a given conductor and order.
    """
    d = 53599
    g = 6

    print(f"Finding the number of primitive Dirichlet characters of conductor d = {d} and order g = {g}.")
    
    # Step 1: Factorize the conductor d.
    prime_factors = get_prime_factorization(d)
    primes = list(prime_factors.keys())
    print(f"Step 1: The prime factorization of d = {d} is {d} = {' * '.join(map(str, sorted(primes)))}.")
    
    # Step 2: Analyze the structure and conditions.
    # For a square-free conductor, primitivity means all components are non-trivial (order > 1).
    # For the final order to be g=6, each component's order must divide 6.
    # Thus, possible orders for each component are {2, 3, 6}.
    print("Step 2: For the character to be primitive and of order 6, each of its components must have an order from the set {2, 3, 6}.")

    # Step 3: Count characters for each possible order using phi(k).
    num_ord_2 = phi(2)
    num_ord_3 = phi(3)
    num_ord_6 = phi(6)
    print("Step 3: The number of characters for each possible order is given by Euler's totient function:")
    print(f"  - Number of characters of order 2: phi(2) = {num_ord_2}")
    print(f"  - Number of characters of order 3: phi(3) = {num_ord_3}")
    print(f"  - Number of characters of order 6: phi(6) = {num_ord_6}")

    # Step 4: Use inclusion-exclusion to find the final count.
    num_components = len(primes)
    choices_per_component = num_ord_2 + num_ord_3 + num_ord_6
    total_combinations = choices_per_component ** num_components
    
    print("Step 4: We use the principle of inclusion-exclusion.")
    print(f"The total number of ways to choose component characters with orders in {{2, 3, 6}} is:")
    print(f"({num_ord_2} + {num_ord_3} + {num_ord_6})^{num_components} = {choices_per_component}^{num_components} = {total_combinations}")

    # Exclude cases where lcm is not 6.
    # Case A: lcm divides 2 (all orders are 2).
    excluded_lcm_div_2 = num_ord_2 ** num_components
    print("We subtract cases where the lcm of orders is not 6.")
    print(f" - Tuples where all orders are 2: {num_ord_2}^{num_components} = {excluded_lcm_div_2}")

    # Case B: lcm divides 3 (all orders are 3).
    excluded_lcm_div_3 = num_ord_3 ** num_components
    print(f" - Tuples where all orders are 3: {num_ord_3}^{num_components} = {excluded_lcm_div_3}")

    # Final calculation
    final_result = total_combinations - excluded_lcm_div_2 - excluded_lcm_div_3
    print("The final number is the total minus the excluded sets.")
    print(f"Final Equation: {total_combinations} - {excluded_lcm_div_2} - {excluded_lcm_div_3} = {final_result}")
    
    print(f"\nThe number of primitive Dirichlet characters is {final_result}.")

solve()
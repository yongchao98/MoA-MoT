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
    """Computes Euler's totient function phi(n)."""
    if n == 1:
        return 1
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
    return int(result)

def solve_dirichlet_problem():
    """
    Calculates the number of primitive Dirichlet characters for a given
    square-free conductor and a specific order.
    """
    d = 53599
    order_k = 6

    print(f"Finding the number of primitive Dirichlet characters of conductor d = {d} and order k = {order_k}.")
    print("-" * 70)

    # Step 1: Factor the conductor d
    prime_factors = get_prime_factorization(d)
    num_prime_factors = len(prime_factors)
    print(f"The conductor d = {d} is square-free with {num_prime_factors} prime factors: {prime_factors}.")

    # Step 2: Set up the counting problem based on character properties
    print("\nA character chi is primitive mod d if it's a product of primitive characters mod p_i.")
    print("A character mod p_i is primitive if its order > 1.")
    print(f"We need lcm(order(chi_1), ..., order(chi_{num_prime_factors})) = {order_k}.")
    print("This means each component's order must be a divisor of 6, and > 1. So, orders can be {2, 3, 6}.")

    # Step 3: Count the number of available characters for each allowed order
    num_chars_order_2 = phi(2)
    num_chars_order_3 = phi(3)
    num_chars_order_6 = phi(6)

    # For each prime factor, count the number of primitive characters with order dividing 6
    num_options_per_prime = num_chars_order_2 + num_chars_order_3 + num_chars_order_6

    print("\nFor each prime factor, we count the number of characters for each allowed order:")
    print(f"Number of characters of order 2 is phi(2) = {num_chars_order_2}")
    print(f"Number of characters of order 3 is phi(3) = {num_chars_order_3}")
    print(f"Number of characters of order 6 is phi(6) = {num_chars_order_6}")
    print(f"Total options per prime factor = {num_chars_order_2} + {num_chars_order_3} + {num_chars_order_6} = {num_options_per_prime}")

    # Step 4: Use inclusion-exclusion to find the final answer
    total_combinations = num_options_per_prime ** num_prime_factors
    
    # Exclude cases where the lcm is not 6
    # Case A: Lcm's prime factors do not include 3 (all orders are 2)
    excluded_A = num_chars_order_2 ** num_prime_factors
    # Case B: Lcm's prime factors do not include 2 (all orders are 3)
    excluded_B = num_chars_order_3 ** num_prime_factors

    final_count = total_combinations - excluded_A - excluded_B

    print("\nUsing inclusion-exclusion to find characters with order exactly 6:")
    print("1. Total combinations of primitive characters with order dividing 6:")
    print(f"   ({num_options_per_prime})^{num_prime_factors} = {total_combinations}")
    print("\n2. Subtract combinations where the order is not 6:")
    print(f"   - Lcm lacks factor 3 (all orders are 2): ({num_chars_order_2})^{num_prime_factors} = {excluded_A}")
    print(f"   - Lcm lacks factor 2 (all orders are 3): ({num_chars_order_3})^{num_prime_factors} = {excluded_B}")
    
    print("\nFinal calculation:")
    print(f"The number of primitive characters is {num_options_per_prime}^{num_prime_factors} - {num_chars_order_2}^{num_prime_factors} - {num_chars_order_3}^{num_prime_factors} = {total_combinations} - {excluded_A} - {excluded_B} = {final_count}")

solve_dirichlet_problem()
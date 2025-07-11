import math

def get_prime_factorization(n):
    """Returns the prime factorization of n as a dictionary."""
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
    """Calculates Euler's totient function."""
    if n == 1:
        return 1
    factors = get_prime_factorization(n)
    result = n
    for p in factors:
        result -= result // p
    return result

def solve():
    """
    Solves the problem of finding the number of primitive Dirichlet characters
    of a given conductor and order.
    """
    d = 53599
    order = 6

    print(f"Finding the number of primitive Dirichlet characters of conductor d = {d} and order k = {order}.")
    
    d_factors = get_prime_factorization(d)
    print(f"Step 1: The prime factorization of the conductor is d = {' * '.join(map(str, sorted(d_factors.keys())))}.")
    
    is_square_free = all(v == 1 for v in d_factors.values())
    if not is_square_free:
        print("The conductor is not square-free. This script only handles square-free conductors.")
        return

    print("The conductor is square-free. A primitive character of conductor d is a product of primitive characters of prime conductors.")
    print("A character mod a prime p is primitive if its order is > 1.")
    print(f"Step 2: The order of the combined character must be lcm(ord(c1), ord(c2), ...) = {order}.")
    print(f"This means component character orders must be divisors of {order} and > 1.")

    k_factors = get_prime_factorization(order)
    prime_factors_of_k = list(k_factors.keys()) # [2, 3] for k=6

    num_p_factors = len(d_factors)

    # Count primitive characters for each prime factor whose order divides 6
    # These orders are 2, 3, 6
    num_ord_2 = phi(2)
    num_ord_3 = phi(3)
    num_ord_6 = phi(6)
    
    num_choices_per_factor = {
        6: num_ord_2 + num_ord_3 + num_ord_6, # Orders dividing 6: {2,3,6}
        3: num_ord_3,                         # Orders dividing 3: {3}
        2: num_ord_2                          # Orders dividing 2: {2}
    }

    print("\nStep 3: Count character choices for each prime factor.")
    print(f"Number of characters of order 2: phi(2) = {num_ord_2}")
    print(f"Number of characters of order 3: phi(3) = {num_ord_3}")
    print(f"Number of characters of order 6: phi(6) = {num_ord_6}")
    
    print("\nStep 4: Use the Principle of Inclusion-Exclusion.")
    # Total combinations where order divides 6
    total_div_6 = num_choices_per_factor[6] ** num_p_factors
    print(f"Number of ways to choose characters where each order divides 6: ({num_ord_2} + {num_ord_3} + {num_ord_6})^{num_p_factors} = {num_choices_per_factor[6]}^{num_p_factors} = {total_div_6}")

    # Combinations to subtract
    # lcm divides 3 (all orders are 3)
    total_div_3 = num_choices_per_factor[3] ** num_p_factors
    print(f"Number of ways where the lcm order divides 3: {num_choices_per_factor[3]}^{num_p_factors} = {total_div_3}")

    # lcm divides 2 (all orders are 2)
    total_div_2 = num_choices_per_factor[2] ** num_p_factors
    print(f"Number of ways where the lcm order divides 2: {num_choices_per_factor[2]}^{num_p_factors} = {total_div_2}")

    result = total_div_6 - total_div_3 - total_div_2
    
    print("\nFinal calculation:")
    print(f"The number of primitive characters of order 6 is ({num_choices_per_factor[6]}^{num_p_factors}) - ({num_choices_per_factor[3]}^{num_p_factors}) - ({num_choices_per_factor[2]}^{num_p_factors}) = {total_div_6} - {total_div_3} - {total_div_2} = {result}")
    
    print("\nThe final answer is:")
    print(result)


solve()
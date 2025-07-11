import math

def get_prime_factorization(n):
    """Returns a list of prime factors of a natural number."""
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
    """Computes Euler's totient function."""
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
    return result

def solve():
    """
    Calculates the number of primitive Dirichlet characters of a given conductor and order.
    """
    d = 53599
    g = 6

    print(f"Finding the number of primitive Dirichlet characters of conductor d = {d} and order g = {g}.\n")

    # Step 1: Factorize the conductor
    prime_factors = get_prime_factorization(d)
    num_factors = len(prime_factors)
    print(f"Step 1: The prime factorization of d = {d} is {prime_factors}.")
    print("Since d is square-free, a character is primitive if and only if it is a product of non-principal characters for each prime factor.\n")
    
    # Step 2: Determine possible orders
    print("Step 2: Determine the possible orders for component characters.")
    print(f"The order of the character is lcm(ord(chi_i), ...) = {g}.")
    print(f"This implies the order of each component character chi_i must be a divisor of {g}.")
    possible_orders = [k for k in range(2, g + 1) if g % k == 0]
    print(f"Since characters must be non-principal (order > 1), the possible orders are {possible_orders}.\n")

    # Step 3: Count characters and use inclusion-exclusion
    print("Step 3: Count the number of valid characters using inclusion-exclusion.")
    
    # Calculate choices for the full set of orders {2, 3, 6}
    choices_per_component = {k: phi(k) for k in possible_orders}
    num_choices_total = sum(choices_per_component.values())
    total_combinations = num_choices_total ** num_factors
    print(f"The number of characters of order 2 is phi(2) = {choices_per_component[2]}.")
    print(f"The number of characters of order 3 is phi(3) = {choices_per_component[3]}.")
    print(f"The number of characters of order 6 is phi(6) = {choices_per_component[6]}.")
    print(f"Total choices per component for an order in {{2, 3, 6}} is {choices_per_component[2]} + {choices_per_component[3]} + {choices_per_component[6]} = {num_choices_total}.")
    print(f"With {num_factors} prime factors, the total combinations where orders are in {{2, 3, 6}} is {num_choices_total}^{num_factors} = {total_combinations}.\n")

    print("Now, we exclude cases where the lcm is not 6.")
    
    # Case A: lcm has no factor of 2. Orders must all be 3.
    num_choices_no_2 = choices_per_component[3]
    count_no_2 = num_choices_no_2 ** num_factors
    print(f"Case A (lcm not divisible by 2): All component orders must be 3.")
    print(f"Number of such combinations = (phi(3))^{num_factors} = {num_choices_no_2}^{num_factors} = {count_no_2}.")

    # Case B: lcm has no factor of 3. Orders must all be 2.
    num_choices_no_3 = choices_per_component[2]
    count_no_3 = num_choices_no_3 ** num_factors
    print(f"Case B (lcm not divisible by 3): All component orders must be 2.")
    print(f"Number of such combinations = (phi(2))^{num_factors} = {num_choices_no_3}^{num_factors} = {count_no_3}.\n")
    
    # Final Calculation
    result = total_combinations - count_no_2 - count_no_3
    print("The final number is the total combinations minus the invalid cases.")
    print(f"Number = {total_combinations} - {count_no_2} - {count_no_3} = {result}")

solve()
<<<608>>>
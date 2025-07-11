import math

def get_prime_factorization(n):
    """
    Returns the prime factorization of n as a dictionary of {prime: exponent}.
    """
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

def solve_dice_problem():
    """
    Calculates the largest possible number of mutually independent events
    for rolling a number of dice with a certain number of sides.
    """
    num_dice = 100
    num_sides = 6

    print(f"Problem: Find the largest number of mutually independent events for rolling {num_dice} {num_sides}-sided dice.")
    print("The sample space size is 6^100 = (2*3)^100 = 2^100 * 3^100.")
    print("This can be decomposed into two independent experiments with sample spaces of size 2^100 and 3^100.")
    print("\nStep 1: Get the prime factorization of the number of sides.")
    prime_factors = get_prime_factorization(num_sides)
    print(f"The prime factorization of {num_sides} is: {prime_factors}")

    print("\nStep 2: For each prime factor p with exponent a, the contribution to the total number of events is a * num_dice.")
    
    contributions = []
    for p, a in prime_factors.items():
        contribution = a * num_dice
        contributions.append(contribution)
        print(f"For prime factor {p} with exponent {a}, the number of independent events is {a} * {num_dice} = {contribution}.")

    total_m = sum(contributions)
    
    print("\nStep 3: The total maximum number of events is the sum of these contributions.")
    
    # Format the final equation string
    equation_parts = [str(c) for c in contributions]
    equation_str = " + ".join(equation_parts)
    
    print("\nFinal Calculation:")
    print(f"{equation_str} = {total_m}")
    
    # The final answer in the required format will be just the number
    # but the user asked for the equation to be printed.
    # The final line of the script will print the final answer.
    print("\nThe largest possible value of m is:")
    print(total_m)


solve_dice_problem()

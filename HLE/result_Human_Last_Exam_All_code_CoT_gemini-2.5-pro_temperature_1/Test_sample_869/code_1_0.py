import math

def calculate_probability():
    """
    Calculates the probability based on the described plan.
    S = Total number of ways to distribute the items.
    F = Number of favorable distributions under the simplifying assumption.
    P = F / S
    """
    
    # Using math.factorial for calculations
    fact5 = math.factorial(5)
    fact25 = math.factorial(25)
    
    # Step 1: Calculate S
    # S = 25! / (5!)^5
    S = fact25 // (fact5**5)
    
    # Step 2: Calculate F
    # Under the simplifying assumption, F is the number of permutations of types.
    # F = 5!
    F = fact5
    
    # Step 3: Calculate P = F / S
    # The probability is F/S. We will print the components of the fraction.
    
    # The problem asks to output each number in the final equation.
    # P = F / S
    print(f"S = {fact25} / ({fact5}^5) = {S}")
    print(f"F = {F}")
    print(f"P = F / S = {F} / {S}")

    # For verification, let's also compute the value of P = (5!)^6 / 25!
    numerator_p = fact5**6
    denominator_p = fact25
    p_value = numerator_p / denominator_p
    
    print(f"The probability can also be expressed as ({fact5}^6) / {fact25}")
    print(f"Value of P is approximately: {p_value}")

calculate_probability()
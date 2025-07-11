import math

def calculate_mistake_bound(n, c):
    """
    Calculates the upper bound on the number of mistakes made by the algorithm.

    Args:
        n (int): The total number of experts.
        c (int): The mistake threshold for removing an expert.

    Returns:
        float: The calculated upper bound.
    """
    # The number of experts n must be at least 1, and c must be at least 1.
    if n < 1 or c < 1:
        print("Number of experts (n) and mistake threshold (c) must be at least 1.")
        return None

    # The formula for the upper bound is c * (n + 1) / 2
    bound = c * (n + 1) / 2
    
    # The problem states the true expert makes strictly fewer than c mistakes,
    # so the number of algorithm mistakes will also be strictly less than this bound.
    # However, the bound itself is a valid upper bound.
    
    print("An upper bound on the number of mistakes is given by the formula: c * (n + 1) / 2")
    # We use math.floor for the result since the number of mistakes is an integer. 
    # The bound M < c*(n+1)/2 means M <= floor(c*(n+1)/2 - epsilon).
    # A simpler integer upper bound is floor(c*(n+1)/2) if c*(n+1) is even, and floor(c*(n+1)/2) if it's odd.
    # An even simpler upper bound is simply the value c*(n+1)/2. The question asks for an upper bound, not the tightest one.
    # Let's present the calculated value.
    
    print(f"For n = {n} and c = {c}:")
    print(f"Bound = {c} * ({n} + 1) / 2 = {bound}")
    
    # A slightly tighter bound is (c-1) + c*(n-1)/2 = c*(n+1)/2 - 1
    tighter_bound = c * (n + 1) / 2 - 1
    print(f"A tighter version of the bound is c*(n+1)/2 - 1 = {tighter_bound}")


# Example values for n and c
n_experts = 11
c_mistakes = 5

calculate_mistake_bound(n_experts, c_mistakes)
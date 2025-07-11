import math

def get_mistake_bound(n, c):
    """
    Calculates the upper bound on the number of mistakes made by the algorithm.

    Args:
        n (int): The total number of experts.
        c (int): The number of mistakes an expert can make before being removed.
    """
    if n <= 1:
        print("The number of experts (n) must be greater than 1.")
        # The true expert is the only expert, so the bound is their mistakes.
        if n == 1:
            print(f"With only one expert (the true one), the bound is c - 1 = {c-1}")
        return

    if c <= 0:
        print("The mistake threshold (c) must be a positive integer.")
        return

    # M_b: Bound for mistakes when the true expert is also wrong.
    # The true expert makes at most c-1 mistakes.
    m_b = c - 1

    # M_g: Bound for mistakes when the true expert is correct.
    # This requires a majority of the n-1 "false" experts to be wrong.
    # We calculate the number of "halving phases" for the false experts.
    
    num_false_experts = n - 1
    phases = 0
    
    f = num_false_experts
    # This loop calculates P(n-1), the number of phases.
    while f > 0:
        phases += 1
        # In each phase, we use the smallest possible majority of the currently
        # active false experts to make c mistakes.
        # The size of this minimal majority is floor(f/2) + 1.
        # After c mistakes, this group is eliminated.
        f = f - (math.floor(f / 2) + 1)
        
    m_g = c * phases
    
    # Total bound is M <= M_g + M_b
    total_bound = m_g + m_b
    
    print(f"For n = {n} experts and a mistake threshold c = {c}:")
    print("The upper bound on the number of mistakes is calculated as:")
    print(f"M <= c * P(n-1) + c - 1")
    print(f"M <= {c} * {phases} + {c} - 1")
    print(f"M <= {m_g} + {m_b}")
    print(f"M <= {total_bound}")

# --- User-configurable values ---
# You can change these values to see the bound for different scenarios.
num_experts = 10
mistake_threshold = 5
# ---

get_mistake_bound(num_experts, mistake_threshold)
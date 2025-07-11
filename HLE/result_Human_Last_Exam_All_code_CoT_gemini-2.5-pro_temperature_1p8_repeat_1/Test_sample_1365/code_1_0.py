import math

def calculate_mistake_bound(n, c):
    """
    Calculates the upper bound on the number of mistakes made by the majority
    voting algorithm given n experts and a mistake limit c.
    
    Args:
        n (int): The total number of experts. Must be >= 2.
        c (int): The number of mistakes an expert must make to be removed.
    
    Returns:
        int: The calculated upper bound for the number of mistakes.
    """
    if not isinstance(n, int) or not isinstance(c, int) or n < 2 or c < 1:
        print("Error: n and c must be integers, with n >= 2 and c >= 1.")
        return None

    # Part 1: Bound for mistakes where the true expert is also wrong (M2)
    # The true expert makes strictly fewer than c mistakes, so M2 <= c - 1.
    m2_bound = c - 1

    # Part 2: Bound for mistakes where the true expert is correct (M1)
    # The minimum number of experts that must be wrong for the majority to err.
    # This happens when no experts have been disqualified yet.
    # All these experts must be "false" experts since the true expert is correct.
    min_wrong_experts_per_mistake = n // 2 + 1
    
    # The total "mistake budget" of all n-1 false experts.
    total_mistake_budget = (n - 1) * c
    
    # M1 is bounded by the total budget divided by the minimum cost per mistake.
    # Since M1 must be an integer, we take the floor.
    m1_bound_float = total_mistake_budget / min_wrong_experts_per_mistake
    m1_bound = math.floor(m1_bound_float)

    # Total bound is the sum of the bounds for M1 and M2.
    total_bound = m1_bound + m2_bound

    # --- Output the results as an equation ---
    print("Derivation of the Upper Bound:")
    print(f"For n = {n} and c = {c}:")
    print(f"\n1. The bound for M2 (algorithm and true expert are wrong) is c - 1.")
    print(f"   M2_bound = {c} - 1 = {m2_bound}")

    print(f"\n2. To find the bound for M1 (algorithm is wrong, true expert is correct):")
    print(f"   a. Total mistake budget from {n-1} false experts = ({n} - 1) * {c} = {total_mistake_budget}")
    print(f"   b. Minimum cost (number of wrong experts) per M1 mistake = floor({n}/2) + 1 = {min_wrong_experts_per_mistake}")
    print(f"   c. M1_bound is floor(budget / cost) = floor({total_mistake_budget} / {min_wrong_experts_per_mistake}) = floor({m1_bound_float:.2f}) = {m1_bound}")

    print(f"\n3. The total upper bound is M1_bound + M2_bound.")
    print(f"   Final Bound = {m1_bound} + {m2_bound} = {total_bound}")
    
    return total_bound

# Example usage with n=11 experts and c=10 mistakes until removal.
calculate_mistake_bound(n=11, c=10)
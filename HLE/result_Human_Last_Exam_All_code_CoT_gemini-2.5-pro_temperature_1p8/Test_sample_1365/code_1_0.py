import math

def calculate_mistake_bound(n, c):
    """
    Calculates and explains the upper bound on mistakes for a majority voting algorithm.

    Args:
      n (int): The total number of experts.
      c (int): The number of mistakes an expert makes before being removed.
    """

    print(f"Calculating the upper bound on algorithm mistakes for n={n} experts and c={c} mistakes for removal.")
    print("-" * 80)
    print("The upper bound is based on the formula: M_total <= M1 + M2\n")

    # Step 1: Explain and calculate the bound for M2
    print("Step 1: Bound M2 (mistakes where the true expert is also wrong)")
    print("A mistake of this type can only happen when the true expert is wrong.")
    print(f"The true expert makes strictly fewer than c mistakes (i.e., at most c-1).")
    m2_bound = c - 1
    print(f"Therefore, M2 <= c - 1")
    print(f"M2 <= {c} - 1 = {m2_bound}\n")

    # Step 2: Explain and calculate the bound for M1
    print("Step 2: Bound M1 (mistakes where the true expert is correct)")
    print("For each of these mistakes, the majority vote was wrong, but the true expert was right.")
    print("This implies that the number of wrong experts was greater than the number of correct experts.")
    print("Since the true expert was correct, all wrong experts must be 'false' experts.")
    print("It can be shown that at least 2 false experts must be wrong for an M1 mistake to occur.")
    
    false_experts = n - 1
    total_false_expert_mistakes = false_experts * c
    print(f"The (n-1) false experts can make at most c mistakes each, for a total budget of (n-1)*c mistakes.")
    print(f"Total mistake budget for false experts = ({n}-1) * {c} = {total_false_expert_mistakes}")
    print(f"Since each M1 mistake costs at least 2 mistakes from this budget, 2 * M1 <= (n-1)*c.")
    
    # M1 must be an integer, so we take the floor of the division.
    m1_bound = math.floor((false_experts * c) / 2)
    print(f"M1 <= floor(((n-1)*c) / 2)")
    print(f"M1 <= floor(({false_experts} * {c}) / 2) = floor({total_false_expert_mistakes / 2}) = {m1_bound}\n")

    # Step 3: Combine the bounds
    print("Step 3: Combine the bounds for M1 and M2")
    total_bound = m1_bound + m2_bound
    print("M_total <= M1 + M2")
    print(f"M_total <= {m1_bound} + {m2_bound}")
    print(f"Final Upper Bound = {total_bound}")
    print("-" * 80)


# --- User-executable code ---
# Example values for n and c
num_experts = 21
mistake_threshold = 10

calculate_mistake_bound(num_experts, mistake_threshold)

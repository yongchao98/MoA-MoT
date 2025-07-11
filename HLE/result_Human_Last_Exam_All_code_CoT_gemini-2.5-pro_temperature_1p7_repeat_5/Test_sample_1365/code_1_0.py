import math

def calculate_mistake_bound(n, c):
    """
    Calculates and prints the upper bound on mistakes for the described majority voting algorithm.

    Args:
        n (int): The total number of experts.
        c (int): The mistake threshold for removing an expert.
    """
    if not isinstance(n, int) or not isinstance(c, int) or n < 1 or c < 1:
        print("Error: n and c must be positive integers.")
        return

    # Part 1: Bound for M_agree
    # M_agree <= c - 1
    m_agree_bound = c - 1

    # Part 2: Bound for M_disagree
    # M_disagree <= floor((n - 1) * c / 2)
    m_disagree_bound = (n - 1) * c // 2 # Using integer division for floor

    # Combine the bounds for the total mistake bound M
    total_bound = m_agree_bound + m_disagree_bound

    # Print the step-by-step derivation with the given numbers
    print(f"For n = {n} experts and mistake threshold c = {c}:")
    print("-------------------------------------------------------")
    print("The upper bound on algorithm mistakes (M) is derived as:")
    print("M <= M_agree_bound + M_disagree_bound\n")

    print("1. Bound for mistakes when the true expert is also wrong (M_agree):")
    print(f"   M_agree <= c - 1")
    print(f"   M_agree <= {c} - 1 = {m_agree_bound}\n")

    print("2. Bound for mistakes when the true expert is right (M_disagree):")
    print(f"   M_disagree <= floor((n - 1) * c / 2)")
    print(f"   M_disagree <= floor(({n} - 1) * {c} / 2) = floor({(n - 1) * c}/2) = {m_disagree_bound}\n")
    
    print("3. Combining these for the total upper bound:")
    print(f"   M <= (c - 1) + floor(((n - 1) * c) / 2)")
    print(f"   M <= {m_agree_bound} + {m_disagree_bound}")
    print(f"   M <= {total_bound}")
    print("-------------------------------------------------------")

# --- Example Execution ---
# You can change these values to see the bound for different scenarios.
num_experts = 11
mistake_threshold = 10

calculate_mistake_bound(num_experts, mistake_threshold)

import math

def calculate_mistake_bound(n, c):
    """
    Calculates the upper bound on the number of mistakes for a majority voting algorithm.

    Args:
        n (int): The total number of experts.
        c (int): The mistake threshold for removing an expert.
    """
    if not isinstance(n, int) or n <= 0:
        print("Error: n (number of experts) must be a positive integer.")
        return
    if not isinstance(c, int) or c <= 0:
        print("Error: c (mistake threshold) must be a positive integer.")
        return

    # The number of mistakes M is bounded by M <= M1 + M2.
    
    # M2 is bounded by the number of mistakes of the true expert (m_true < c).
    # M2 <= c - 1
    m2_bound = c - 1

    # M1 is bounded by half the total mistake capacity of the n-1 other experts.
    # M1 <= floor((n - 1) * c / 2)
    m1_bound = math.floor((n - 1) * c / 2)

    # The total upper bound is the sum of the two bounds.
    total_bound = m1_bound + m2_bound

    # Print the final equation with the numbers plugged in.
    print(f"For n = {n} experts and a mistake threshold c = {c}:")
    print("The upper bound on the number of mistakes (M) is calculated as:")
    print(f"M <= floor(({n} - 1) * {c} / 2) + ({c} - 1)")
    print(f"M <= floor({(n - 1) * c / 2}) + {c - 1}")
    print(f"M <= {m1_bound} + {m2_bound}")
    print(f"M <= {total_bound}")

# --- Example Parameters ---
# You can change these values to see the bound for different scenarios.
num_experts = 21
mistake_threshold = 10
# --- End of Parameters ---

calculate_mistake_bound(num_experts, mistake_threshold)
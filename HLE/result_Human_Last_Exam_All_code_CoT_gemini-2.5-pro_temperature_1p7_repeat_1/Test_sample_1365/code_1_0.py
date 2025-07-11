import math

def calculate_mistake_bound(n, c):
    """
    Calculates the upper bound on the number of mistakes for the experts problem variant.

    Args:
        n (int): The total number of experts. Must be >= 1.
        c (int): The number of mistakes an expert can make before being removed. Must be >= 1.
    """
    if n < 1 or c < 1:
        print("Number of experts (n) and mistake threshold (c) must be at least 1.")
        return

    # M_B is the bound for mistakes when the true expert is also wrong.
    m_b_bound = c - 1

    # M_A is the bound for mistakes when the true expert is correct.
    # This is analogous to c parallel runs of the Halving algorithm on the n-1 false experts.
    # If n=1 or n=2, log2(n-1) is <= 0, so the floor will be <= 0.
    # A M_A mistake is not possible if there is only one false expert, so M_A should be 0.
    # We handle the n=1 case to avoid math.log2(0).
    if n <= 1:
        log_term_val = 0 # No false experts to cause an M_A mistake
    else:
        # For n > 1, n-1 is at least 1, so log2 is non-negative.
        log_term_val = math.floor(math.log2(n - 1))
    
    m_a_bound = c * log_term_val
    
    # Total mistake bound is the sum of the two bounds.
    total_bound = m_a_bound + m_b_bound

    print("The upper bound on the number of mistakes is derived from M <= M_A + M_B")
    print(f"Bound on M_B (true expert is wrong): {m_b_bound}")
    print(f"Bound on M_A (true expert is right): c * floor(log2(n-1)) = {c} * {log_term_val} = {m_a_bound}")
    print("Final Upper Bound M <= (c-1) + c * floor(log2(n-1))")
    print(f"M <= ({c} - 1) + {c} * floor(log2({n} - 1))")
    print(f"M <= {m_b_bound} + {m_a_bound}")
    print(f"M <= {total_bound}")

# Example usage with n=101 experts and c=20 mistake threshold.
# You can change these values to see the bound for different scenarios.
n_experts = 101
mistake_threshold = 20
calculate_mistake_bound(n_experts, mistake_threshold)
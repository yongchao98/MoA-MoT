import math

def calculate_mistake_bound(n, c):
    """
    Calculates the upper bound on the number of mistakes for a variant of the experts problem.

    Args:
        n (int): The total number of experts.
        c (int): The number of mistakes after which an expert is removed.

    Returns:
        int: The upper bound on the number of mistakes made by the algorithm.
    """
    if not isinstance(n, int) or n <= 0:
        raise ValueError("Number of experts 'n' must be a positive integer.")
    if not isinstance(c, int) or c <= 0:
        raise ValueError("Mistake threshold 'c' must be a positive integer.")

    # Bound for mistakes where the true expert is also wrong.
    # M_B <= c - 1
    m_b_bound = c - 1

    # Bound for mistakes where the true expert is correct.
    # M_A <= floor(c * (n - 1) / 2)
    m_a_bound = math.floor(c * (n - 1) / 2)

    # Total mistake bound M = M_A + M_B
    total_bound = m_b_bound + m_a_bound

    # Output the equation with the numbers plugged in
    print("Deriving the upper bound for n={} and c={}:".format(n, c))
    print("Bound = (c - 1) + floor(c * (n - 1) / 2)")
    print("Bound = ({} - 1) + floor({} * ({} - 1) / 2)".format(c, c, n))
    print("Bound = {} + floor({} / 2)".format(m_b_bound, c * (n - 1)))
    print("Bound = {} + {}".format(m_b_bound, m_a_bound))
    print("Upper Bound on Mistakes = {}".format(total_bound))

    return total_bound

# --- Example Usage ---
# You can change these values to see the bound for different scenarios.
num_experts = 11
mistake_threshold = 10

calculate_mistake_bound(num_experts, mistake_threshold)
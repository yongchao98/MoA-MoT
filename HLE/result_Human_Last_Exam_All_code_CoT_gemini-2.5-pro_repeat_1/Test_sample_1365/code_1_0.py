import math

def calculate_mistake_bound(n, c):
    """
    Calculates the upper bound on the number of mistakes for a majority voting algorithm.

    Args:
        n (int): The total number of experts.
        c (int): The number of mistakes after which an expert is removed.
    """
    if not isinstance(n, int) or not isinstance(c, int) or n < 1 or c < 1:
        print("Error: 'n' and 'c' must be positive integers.")
        return

    # Bound for M_W (mistakes when the true expert is wrong).
    # This is at most the number of mistakes the true expert can make.
    # m_true < c  =>  m_true <= c - 1
    bound_m_w = c - 1

    # Bound for M_R (mistakes when the true expert is right).
    # Each such mistake requires at least 2 mistakes from the non-true experts.
    # The total mistake budget for non-true experts is (n-1)*c.
    # So, 2 * M_R <= (n-1)*c. We use integer division for the floor.
    bound_m_r = (n - 1) * c // 2

    # The total upper bound is the sum of the two parts.
    upper_bound = bound_m_w + bound_m_r

    # Print the explanation and the final equation.
    print(f"For n = {n} experts and mistake cutoff c = {c}:")
    print("The upper bound on algorithm mistakes (M) is calculated as M <= M_W + M_R.")
    print(f"M_W <= c - 1 = {c} - 1 = {bound_m_w}")
    print(f"M_R <= floor((n - 1) * c / 2) = floor(({n} - 1) * {c} / 2) = {bound_m_r}")
    print("--------------------------------------------------")
    print(f"Total Mistakes M <= {bound_m_w} + {bound_m_r} = {upper_bound}")

# --- Example Usage ---
# You can change these values to see the bound for different scenarios.
num_experts = 11
mistake_cutoff = 10

calculate_mistake_bound(num_experts, mistake_cutoff)

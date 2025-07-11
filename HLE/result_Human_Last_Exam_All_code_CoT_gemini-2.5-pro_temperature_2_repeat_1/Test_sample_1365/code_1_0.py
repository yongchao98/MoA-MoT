import math

def calculate_mistake_bound(n, c):
    """
    Calculates the upper bound on the number of mistakes made by the algorithm.

    Args:
      n: The total number of experts.
      c: The number of mistakes an expert must make to be removed.
    """
    if n <= 1:
        # If there's only one expert, they must be the true expert.
        # The algorithm will just follow them and make the same number of mistakes.
        # This is at most c-1.
        bound = c - 1
        print(f"With n=1, the bound is simply c-1 = {c-1}")
        print(f"<<<answer {c-1}>>>")
        return

    # M_b: Mistakes when the true expert is wrong. This is at most c-1.
    m_b = c - 1

    # M_g: Mistakes when the true expert is correct.
    # The number of "halving epochs" for the n-1 non-true experts.
    log_term = math.floor(math.log2(n - 1))
    
    # Each epoch takes at most c mistakes.
    m_g = c * log_term

    # Total mistake bound M = M_b + M_g
    total_bound = m_b + m_g
    
    print(f"An upper bound on the number of mistakes is given by the formula:")
    print(f"M <= (c - 1) + c * floor(log2(n - 1))")
    print(f"M <= ({c} - 1) + {c} * floor(log2({n} - 1))")
    print(f"M <= {m_b} + {c} * {log_term}")
    print(f"M <= {m_b} + {m_g}")
    print(f"M <= {total_bound}")
    print(f"\nFinal calculated upper bound:")
    print(total_bound)
    print(f"<<<answer {total_bound}>>>")


# --- User Input ---
# You can change these values to see the bound for different scenarios.
n_experts = 10  # Total number of experts
c_mistakes = 5  # Mistake threshold for removal

calculate_mistake_bound(n_experts, c_mistakes)
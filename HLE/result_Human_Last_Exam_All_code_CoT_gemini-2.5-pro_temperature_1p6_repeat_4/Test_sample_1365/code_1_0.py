def calculate_mistake_bound(n, c):
    """
    Calculates the upper bound on the number of mistakes made by the algorithm.

    Args:
        n (int): The total number of experts.
        c (int): The number of mistakes an expert can make before being removed.
    """
    if n <= 1 or c <= 0:
        print("n must be greater than 1 and c must be greater than 0.")
        return

    # The number of mistakes of the true expert, m_star, is strictly less than c.
    # So, the maximum value for m_star is c - 1.
    max_m_star = c - 1

    # The number of non-true experts is n - 1.
    # The total mistake budget for these experts is (n - 1) * c.
    # We found that D_R (number of algorithm mistakes when the true expert is right)
    # is less than or equal to this budget.
    dr_bound = (n - 1) * c

    # The total number of algorithm mistakes M is bounded by D_R + m_star.
    upper_bound = dr_bound + max_m_star

    print(f"Given n = {n} experts and a mistake limit of c = {c}:")
    print("\nAn upper bound on the number of algorithm mistakes M is derived from:")
    print("M <= (Number of non-true experts) * c + (mistakes of true expert)")
    print(f"M <= (n-1)*c + m_star")
    print(f"Since m_star < c, the maximum integer value for m_star is c-1 = {max_m_star}.")
    print("\nSubstituting the maximum value for m_star gives the upper bound:")
    # The final equation as requested by the user prompt
    print(f"M <= ({n}-1)*{c} + ({c}-1)")
    print(f"M <= {n-1}*{c} + {max_m_star}")
    print(f"M <= {dr_bound} + {max_m_star}")
    print(f"M <= {upper_bound}")


# Example values for n and c
n_experts = 10
c_mistakes = 50

calculate_mistake_bound(n_experts, c_mistakes)
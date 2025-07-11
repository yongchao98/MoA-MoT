def calculate_mistake_bound(n, c):
    """
    Calculates the upper bound on the number of mistakes for a variant of the experts problem.

    Args:
        n (int): The total number of experts.
        c (int): The mistake threshold for removing an expert.
    """
    if not isinstance(n, int) or n <= 0:
        print("Error: n must be a positive integer.")
        return
    if not isinstance(c, int) or c <= 0:
        print("Error: c must be a positive integer.")
        return

    # The formula for the upper bound M is: M <= (c - 1) + (n - 1) * c / 2
    
    # Bound for mistakes when the true expert is also wrong.
    term1 = c - 1
    
    # Bound for mistakes when the true expert is right.
    term2_numerator = (n - 1) * c
    term2_denominator = 2
    term2 = term2_numerator / term2_denominator
    
    # Total upper bound
    upper_bound = term1 + term2

    print("An upper bound on the number of mistakes (M) is given by the formula:")
    print("M <= (c - 1) + (n - 1) * c / 2")
    print(f"\nFor n = {n} and c = {c}, the calculation is:")
    print(f"M <= ({c} - 1) + ({n} - 1) * {c} / {term2_denominator}")
    print(f"M <= {term1} + {term2_numerator} / {term2_denominator}")
    print(f"M <= {term1} + {term2}")
    print(f"M <= {upper_bound}")


# --- Example Usage ---
# Total number of experts
n_experts = 21
# Mistake threshold for removing an expert
mistake_cutoff = 10

calculate_mistake_bound(n_experts, mistake_cutoff)

import math

def calculate_mistake_bound():
    """
    Calculates and prints an upper bound on the number of mistakes made by a
    majority voting algorithm in a specific variant of the experts problem.
    """
    # --- Configuration ---
    # n: total number of experts (must be > 1)
    n = 20
    # c: mistake limit after which an expert is removed (must be >= 1)
    # The true expert makes strictly fewer than c mistakes.
    c = 10
    # ---------------------

    if n <= 1:
        print("Error: The number of experts (n) must be greater than 1.")
        return
    if c < 1:
        print("Error: The mistake limit (c) must be 1 or greater.")
        return

    # The upper bound on the number of mistakes 'M' is given by the formula:
    # M <= c * log2(n-1) + (c-1)
    # This is derived by splitting the algorithm's mistakes into two cases:
    # 1. When the true expert is wrong (at most c-1 times).
    # 2. When the true expert is correct. In this case, a majority of the other
    #    n-1 experts must be wrong. The number of such events is bounded by
    #    c * log2(n-1), a standard result for this type of "sleeping experts" problem.

    # Calculate the components of the bound
    log_val = math.log2(n - 1)
    c_minus_1 = c - 1
    bound = c * log_val + c_minus_1

    # Print the formula and the step-by-step calculation
    print(f"An upper bound on the number of mistakes (M) is given by the formula:")
    print(f"M <= c * log2(n - 1) + (c - 1)")
    print(f"\nFor n = {n} and c = {c}, the calculation is as follows:")
    print(f"M <= {c} * log2({n - 1}) + ({c} - 1)")
    print(f"M <= {c} * {log_val:.4f} + {c_minus_1}")
    print(f"M <= {c * log_val:.4f} + {c_minus_1}")
    print(f"M <= {bound:.4f}")

calculate_mistake_bound()
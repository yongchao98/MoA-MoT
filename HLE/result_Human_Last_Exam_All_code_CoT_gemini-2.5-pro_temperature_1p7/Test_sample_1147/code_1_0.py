import math

def calculate_loo_coverage_guarantee(n, alpha):
    """
    Calculates and explains the lowest possible coverage probability for the
    leave-one-out (LOO) conformal prediction scheme.

    Args:
      n (int): The number of training data points.
      alpha (float): The desired miscoverage rate (e.g., 0.1 for 90% confidence).
    """

    print(f"Given n = {n} training points and a desired miscoverage rate alpha = {alpha}.")
    print("The final coverage probability P is determined by the formula: P = k / (n+1)")
    print("where k = ceil((1-alpha) * (n+1)).\n")

    # --- Step-by-step calculation ---
    print("Let's calculate the values for the equation:")

    # 1. Calculate n+1
    n_plus_1 = n + 1
    print(f"1. The total number of LOO scores (including the test point) is n+1.")
    print(f"   n + 1 = {n} + 1 = {n_plus_1}")

    # 2. Calculate 1-alpha
    one_minus_alpha = 1 - alpha
    print(f"2. The target coverage level is 1-alpha.")
    print(f"   1 - alpha = 1 - {alpha} = {one_minus_alpha}")

    # 3. Calculate k
    k_float = one_minus_alpha * n_plus_1
    k = math.ceil(k_float)
    print(f"3. The quantile rank 'k' is the ceiling of their product.")
    print(f"   k = ceil({one_minus_alpha} * {n_plus_1}) = ceil({k_float:.2f}) = {k}")

    # 4. Calculate the final probability
    lowest_coverage = k / n_plus_1
    print(f"4. The lowest possible coverage is k / (n+1).")
    print(f"   P = {k} / {n_plus_1} = {lowest_coverage:.4f}\n")

    print(f"Therefore, for any data distribution and any algorithm, the coverage probability")
    print(f"is guaranteed to be at least {lowest_coverage:.4f} (or {lowest_coverage:.2%}).")

# --- Execute with example values ---
# You can change these values to see the result for different scenarios.
example_n = 18
example_alpha = 0.1

calculate_loo_coverage_guarantee(example_n, example_alpha)
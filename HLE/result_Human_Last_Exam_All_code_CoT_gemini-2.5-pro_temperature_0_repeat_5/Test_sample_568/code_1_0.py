import math

def calculate_mixing_repetitions(epsilon, delta):
    """
    Calculates the number of repetitions needed for the liquids in two cups to become "the same".

    Args:
        epsilon (float): The maximum allowed difference in concentration (must be between 0 and 1).
        delta (float): The fraction of liquid exchanged in each step (must be between 0 and 1).
    """
    # Validate inputs
    if not (0 < epsilon < 1):
        print("Error: Epsilon (ε) must be in the interval (0, 1).")
        return
    if not (0 < delta < 1):
        print("Error: Delta (δ) must be in the interval (0, 1).")
        return

    # Handle the special case where delta = 0.5
    if delta == 0.5:
        n = 1
        print("For delta = 0.5, perfect mixing is achieved in a single step.")
        print(f"The number of repetitions needed is {n}.")
        return

    # Calculate the terms of the formula
    log_epsilon = math.log(epsilon)
    base = abs(1 - 2 * delta)
    log_base = math.log(base)
    
    # Calculate n
    ratio = log_epsilon / log_base
    n = math.ceil(ratio)

    # Print the final equation with all numbers substituted
    print("The number of repetitions 'n' is calculated using the formula:")
    print("n = ceil(ln(ε) / ln(|1 - 2δ|))\n")
    print("Substituting the given values:")
    final_equation = f"n = ceil(ln({epsilon}) / ln(|1 - 2*{delta}|))"
    final_calculation = f"n = ceil({log_epsilon:.4f} / {log_base:.4f})"
    final_result = f"n = ceil({ratio:.4f}) = {n}"
    
    print(final_equation)
    print(final_calculation)
    print(final_result)


# --- Parameters ---
# You can change these values to see the result for different parameters.
epsilon = 0.01  # Example value for ε
delta = 0.1     # Example value for δ

# --- Execution ---
calculate_mixing_repetitions(epsilon, delta)
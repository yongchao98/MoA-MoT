import math

def calculate_mixing_steps(epsilon, delta):
    """
    Calculates the number of steps to mix two liquids based on epsilon and delta.

    Args:
        epsilon (float): The maximum allowed difference in concentration (0 < epsilon < 1).
        delta (float): The fraction of liquid moved in each step (0 < delta < 1).
    """
    if not (0 < epsilon < 1 and 0 < delta < 1):
        print("Error: Epsilon and delta must be between 0 and 1.")
        return

    # Handle edge cases where mixing never completes.
    if delta == 0 or delta == 1:
        print(f"For delta = {delta}, the liquids will never reach the desired mixed state.")
        print("The number of steps is effectively infinite.")
        return

    # Handle the fastest mixing case.
    if delta == 0.5:
        n = 1
        print("For delta = 0.5, mixing is completed in a single step.")
        print(f"The number of steps required is {n}.")
        return

    print(f"Calculating the number of steps for epsilon = {epsilon} and delta = {delta}:")
    print("\nThe formula is n = ceil(log(epsilon) / log(|1 - 2*delta|))")
    print("---------------------------------------------------------")

    # Step 1: Calculate the base of the exponent
    base = abs(1 - 2 * delta)
    print(f"1. The base of the exponent is |1 - 2*({delta})| = {base}")

    # Step 2: Calculate the logarithms
    log_epsilon = math.log(epsilon)
    log_base = math.log(base)
    print(f"2. log(epsilon) = log({epsilon}) = {log_epsilon:.4f}")
    print(f"   log(base) = log({base}) = {log_base:.4f}")

    # Step 3: Calculate the ratio
    ratio = log_epsilon / log_base
    print(f"3. The ratio is {log_epsilon:.4f} / {log_base:.4f} = {ratio:.4f}")

    # Step 4: Find the smallest integer n
    n = math.ceil(ratio)
    print(f"4. The smallest integer n >= {ratio:.4f} is ceil({ratio:.4f}) = {n}")
    print("---------------------------------------------------------")
    print(f"Result: The procedure must be repeated {n} times.")


# --- Example Usage ---
# You can change these values to see the result for different parameters.
epsilon_value = 0.05  # We want the liquids to be 95% the same
delta_value = 0.1   # We exchange 10% of the liquid at each step

calculate_mixing_steps(epsilon_value, delta_value)
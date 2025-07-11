import math

def calculate_mixing_steps(epsilon, delta):
    """
    Calculates the number of steps to mix two liquids until they are "the same".

    Args:
        epsilon (float): The maximum allowed difference in concentration (0 < epsilon < 1).
        delta (float): The fraction of liquid swapped in each step (0 < delta < 1).
    """
    print(f"Calculating for epsilon = {epsilon} and delta = {delta}")
    print("-----------------------------------------")

    if not (0 < epsilon < 1):
        print("Error: Epsilon must be between 0 and 1.")
        return
    if not (0 < delta < 1):
        print("Error: Delta must be between 0 and 1.")
        return

    # Handle the special case where mixing is perfect in one step
    if delta == 0.5:
        n = 1
        print("When delta is 0.5, the liquids mix perfectly in a single step.")
        print(f"Final answer n = {n}")
        return

    # The formula is n = ceil( log(epsilon) / log(|1 - 2*delta|) )
    # Output the numbers in the final equation as requested
    
    val = abs(1 - 2 * delta)
    log_epsilon = math.log(epsilon)
    log_abs_val = math.log(val)
    ratio = log_epsilon / log_abs_val
    n = math.ceil(ratio)

    print(f"The number of steps 'n' is the smallest integer satisfying |1 - 2*delta|^n <= epsilon")
    print(f"We need to solve for n in: |1 - 2*{delta}|^n <= {epsilon}")
    print(f"This is equivalent to: n >= log({epsilon}) / log(|1 - 2*{delta}|)")
    print(f"n >= log({epsilon}) / log({val})")
    print(f"n >= {log_epsilon} / {log_abs_val}")
    print(f"n >= {ratio}")
    print(f"Since n must be an integer, we take the ceiling.")
    print(f"n = ceil({ratio})")
    print(f"Final answer n = {n}")


# Example usage:
epsilon_val = 0.01
delta_val = 0.1
calculate_mixing_steps(epsilon_val, delta_val)
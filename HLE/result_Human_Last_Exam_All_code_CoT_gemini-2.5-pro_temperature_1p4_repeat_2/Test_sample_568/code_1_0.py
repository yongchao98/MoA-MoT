import math

def calculate_mixing_steps(epsilon, delta):
    """
    Calculates the number of steps required for the liquids in two cups to be considered "the same".

    Args:
        epsilon (float): The maximum allowed difference in concentration (0 < epsilon < 1).
        delta (float): The fraction of liquid exchanged in each step (0 < delta < 1).

    Returns:
        None. Prints the detailed calculation.
    """
    # Print problem parameters
    print(f"Given parameters: epsilon = {epsilon}, delta = {delta}")
    print("-" * 30)

    # The final formula for the number of steps t is:
    # t = ceil(log(epsilon) / log(|1 - 2*delta|))

    # Handle the special case where delta = 0.5
    if delta == 0.5:
        t = 1
        print("Special case: delta = 0.5")
        print("The liquids mix perfectly in one step.")
        print("Final Answer, t = 1")
        return

    # Handle invalid delta values that would lead to infinite steps
    if delta <= 0 or delta >= 1:
        print("Invalid delta: delta must be in the range (0, 1).")
        print("For delta=0 or delta=1, the liquids never become the same.")
        return

    # General case calculation
    print("The number of repetitions t is given by the formula:")
    print("t = ceil(log(epsilon) / log(|1 - 2*delta|))")
    print("\nCalculation steps:")
    
    # Numerator of the argument to ceil
    numerator = math.log(epsilon)
    # Denominator of the argument to ceil
    base = abs(1 - 2 * delta)
    denominator = math.log(base)
    
    # The ratio
    ratio = numerator / denominator
    
    # The final number of steps
    t = math.ceil(ratio)

    print(f"1. Calculate the argument of ceil: log({epsilon}) / log(|1 - 2*{delta}|)")
    print(f"   = log({epsilon}) / log(|{1 - 2 * delta}|)")
    print(f"   = {numerator:.4f} / log({base:.4f})")
    print(f"   = {numerator:.4f} / {denominator:.4f}")
    print(f"   = {ratio:.4f}")
    print(f"\n2. Apply the ceiling function: ceil({ratio:.4f})")
    print(f"   = {t}")
    print("-" * 30)
    print(f"Final Answer, t = {t}")


# Example usage with some parameters
# You can change these values to see the result for different epsilon and delta
epsilon_param = 0.01
delta_param = 0.1
calculate_mixing_steps(epsilon_param, delta_param)
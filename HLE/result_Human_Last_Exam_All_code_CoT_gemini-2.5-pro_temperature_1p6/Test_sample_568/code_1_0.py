import math

def calculate_mixing_steps(epsilon, delta):
    """
    Calculates the number of steps required for the liquids in two cups to mix
    until their compositions differ by at most epsilon.

    Args:
        epsilon (float): The tolerance for difference, must be in (0, 1).
        delta (float): The fraction of liquid moved, must be in (0, 1).
    """

    print(f"Calculating steps for epsilon = {epsilon} and delta = {delta}")
    print("--------------------------------------------------")

    # Input validation
    if not (0 < epsilon < 1 and 0 < delta < 1):
        print("Error: Epsilon and delta must be in the interval (0, 1).")
        return

    # Handle the special case where delta = 0.5
    if delta == 0.5:
        n = 1
        print("Special case: delta = 0.5")
        print("The liquids mix perfectly in one step.")
        print("Final Answer n = 1")
    else:
        # The number of steps 'n' is given by the formula:
        # n = ceil( log(epsilon) / log(|1 - 2*delta|) )
        print("The number of steps 'n' is the ceiling of the following equation:")
        print("n = log(epsilon) / log(|1 - 2*delta|)\n")

        # Calculate each part of the formula and print it
        log_epsilon = math.log(epsilon)
        
        base_of_power = 1 - 2 * delta
        abs_base = abs(base_of_power)
        log_abs_base = math.log(abs_base)
        
        ratio = log_epsilon / log_abs_base
        n = math.ceil(ratio)

        print("--- Calculation Breakdown ---")
        print(f"The equation is: ceil( log({epsilon}) / log(|1 - 2*{delta}|) )")
        print(f"1. Numerator: log({epsilon}) = {log_epsilon:.4f}")
        print(f"2. Denominator: log(|{base_of_power:.2f}|) = log({abs_base:.2f}) = {log_abs_base:.4f}")
        print(f"3. Ratio: {log_epsilon:.4f} / {log_abs_base:.4f} = {ratio:.4f}")
        print(f"4. Ceiling (smallest integer n): ceil({ratio:.4f}) = {n}")
        print("--------------------------------------------------")
        print(f"Final Answer: The procedure must be repeated {n} times.")

# --- Example Usage ---
# You can change these values to see the result for different parameters.
epsilon_param = 0.01
delta_param = 0.1

calculate_mixing_steps(epsilon_param, delta_param)
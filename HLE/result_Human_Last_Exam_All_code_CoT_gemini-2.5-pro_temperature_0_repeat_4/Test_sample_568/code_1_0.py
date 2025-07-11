import math

def solve_mixing_problem(epsilon, delta):
    """
    Calculates the number of repetitions needed for the liquids in two cups to be "the same".

    Args:
        epsilon (float): The maximum allowed difference fraction, between 0 and 1.
        delta (float): The fraction of liquid exchanged in each step, between 0 and 1.
    """
    # --- Input Validation ---
    if not (0 < epsilon < 1):
        print("Error: Epsilon (ε) must be between 0 and 1.")
        return
    if not (0 < delta < 1):
        print("Error: Delta (δ) must be between 0 and 1.")
        return

    print(f"Given parameters: ε = {epsilon}, δ = {delta}\n")

    # --- Explanation ---
    print("Let r_A(t) be the fraction of red liquid in cup A after t steps.")
    print("The state of the system evolves according to the recurrence relation:")
    print(f"r_A(t+1) = (1 - 2*δ) * r_A(t) + δ")
    print("The closed-form solution for the fraction of red liquid in cup A is:")
    print(f"r_A(t) = 1/2 + 1/2 * (1 - 2*δ)^t\n")
    
    print("The condition for the liquids to be 'the same' is that the difference")
    print("in their compositions is at most ε. This can be expressed as:")
    print(f"|1 - 2*δ|^t ≤ ε\n")

    print("To find the number of repetitions 't', we solve this inequality for the smallest integer t.")
    print("Taking the logarithm of both sides gives the formula:")
    print("t ≥ log(ε) / log(|1 - 2*δ|)")
    print("So, t = ceil(log(ε) / log(|1 - 2*δ|))\n")

    # --- Calculation ---
    print("Plugging in the values:")
    # Handle the special case where delta = 0.5
    if delta == 0.5:
        print("When δ = 0.5, perfect mixing occurs in a single step.")
        t = 1
        print(f"The required number of repetitions is: {t}")
        return

    val_1_minus_2d = 1 - 2 * delta
    abs_val = abs(val_1_minus_2d)
    
    log_numerator = math.log(epsilon)
    log_denominator = math.log(abs_val)
    
    ratio = log_numerator / log_denominator
    result = math.ceil(ratio)

    print(f"t = ceil(log({epsilon}) / log(|1 - 2*{delta}|))")
    print(f"t = ceil(log({epsilon}) / log(|{val_1_minus_2d}|))")
    print(f"t = ceil(log({epsilon}) / log({abs_val}))")
    print(f"t = ceil({log_numerator} / {log_denominator})")
    print(f"t = ceil({ratio})")
    print(f"\nThe required number of repetitions is: {result}")


# --- Example Usage ---
# You can change these values to see the result for different parameters.
epsilon_param = 0.01
delta_param = 0.1

solve_mixing_problem(epsilon_param, delta_param)
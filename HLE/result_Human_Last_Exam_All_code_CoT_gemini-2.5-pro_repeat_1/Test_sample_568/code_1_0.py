import math

def solve_mixing_problem(epsilon, delta):
    """
    Calculates the number of repetitions (n) needed for the cups to be "the same".

    Args:
        epsilon (float): The maximum allowed difference in concentration, in (0, 1).
        delta (float): The fraction of liquid moved in each step, in (0, 1).
    """
    print(f"Given parameters:\n  epsilon = {epsilon}\n  delta = {delta}\n")

    if not (0 < epsilon < 1 and 0 < delta < 1):
        print("Error: Epsilon and delta must be in the range (0, 1).")
        return

    # Handle the special case where delta = 0.5
    if delta == 0.5:
        n = 1
        print("Special case: delta = 0.5. Mixing is perfect in one step.")
        print(f"Result: n = {n}")
        return

    # The general formula
    # n >= log(epsilon) / log(|1 - 2*delta|)
    log_epsilon = math.log(epsilon)
    base = abs(1 - 2 * delta)
    log_base = math.log(base)
    
    # Calculate n
    n_float = log_epsilon / log_base
    n = math.ceil(n_float)
    
    # Output the calculation step-by-step
    print("The required number of steps 'n' is given by the formula:")
    print("n = ceil(log(epsilon) / log(|1 - 2*delta|))\n")
    
    print("Substituting the given values:")
    print(f"n = ceil(log({epsilon}) / log(|1 - 2*{delta}|))")
    print(f"n = ceil({log_epsilon:.4f} / log({base:.4f}))")
    print(f"n = ceil({log_epsilon:.4f} / {log_base:.4f})")
    print(f"n = ceil({n_float:.4f})")
    print(f"\nResult: n = {n}")


# --- Example Usage ---
# You can change these values to test other scenarios
epsilon_val = 0.01
delta_val = 0.1

solve_mixing_problem(epsilon_val, delta_val)
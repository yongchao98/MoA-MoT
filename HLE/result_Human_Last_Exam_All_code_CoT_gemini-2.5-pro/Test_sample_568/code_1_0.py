import math

def solve_mixing_problem(epsilon, delta):
    """
    Calculates the number of steps required for the liquids to be considered "the same".

    Args:
        epsilon (float): The tolerance for the difference in concentrations, between 0 and 1.
        delta (float): The fraction of liquid exchanged in each step, between 0 and 1.
    """
    print(f"Goal: Find the smallest integer n such that the liquids in the cups are 'the same'.")
    print("This is interpreted as the point where the difference in the concentration of red liquid")
    print("between the two cups is at most epsilon.")
    print(f"\nParameters: epsilon = {epsilon}, delta = {delta}")
    
    print("\nThe derived inequality to solve for n is: |1 - 2*delta|^n <= epsilon")
    print("This leads to the formula: n = ceil(log(epsilon) / log(|1 - 2*delta|))")
    print("-" * 20)

    if not (0 < epsilon < 1 and 0 < delta < 1):
        print("Error: Epsilon and delta must be in the interval (0, 1).")
        return

    if delta == 0.5:
        n = 1
        print("Special case: delta = 0.5")
        print("Mixing is perfectly balanced and completes in a single step.")
        print(f"Final answer: n = {n}")
    else:
        # Calculate the terms of the formula
        base = abs(1 - 2 * delta)
        log_epsilon = math.log(epsilon)
        log_base = math.log(base)
        
        # Perform the division and ceiling operation
        n_float = log_epsilon / log_base
        n = math.ceil(n_float)

        # Print the step-by-step calculation
        print("Substituting the values into the formula:")
        print(f"n = ceil(log({epsilon}) / log(|1 - 2*{delta}|))")
        print(f"n = ceil({log_epsilon:.4f} / log({base:.4f}))")
        print(f"n = ceil({log_epsilon:.4f} / {log_base:.4f})")
        print(f"n = ceil({n_float:.4f})")
        print(f"Final answer: n = {int(n)}")

# --- Example Usage ---
# You can change these values to see the result for different parameters.
epsilon_param = 0.01
delta_param = 0.25

solve_mixing_problem(epsilon_param, delta_param)

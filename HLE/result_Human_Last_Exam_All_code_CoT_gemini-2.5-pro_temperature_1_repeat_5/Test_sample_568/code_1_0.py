import math

def solve_mixing_problem():
    """
    Calculates the number of steps required to mix two liquids until they are
    considered "the same" based on the given parameters epsilon and delta.
    """
    # Fix parameters for the calculation.
    # Epsilon is the maximum allowed difference in composition.
    # Delta is the fraction of liquid exchanged in each step.
    epsilon = 0.01
    delta = 0.1

    print(f"Finding the number of steps 'n' for given epsilon={epsilon} and delta={delta}.")
    print("-" * 20)

    # Handle edge cases for delta.
    if not 0 < delta < 1:
        print("Error: delta must be a value between 0 and 1.")
        return
    if not 0 < epsilon < 1:
        print("Error: epsilon must be a value between 0 and 1.")
        return

    # Handle the special case where delta = 0.5, which leads to perfect mixing in one step.
    if delta == 0.5:
        n = 1
        print("For delta = 0.5, perfect mixing is achieved in a single step.")
        print("n = 1")
        return

    # The general formula is derived from the inequality: |1 - 2*delta|^n <= epsilon
    # This leads to n >= log(epsilon) / log(|1 - 2*delta|)
    # Since n must be an integer, we take the ceiling of the result.
    
    # Calculation steps
    base = abs(1 - 2 * delta)
    log_epsilon = math.log(epsilon)
    log_base = math.log(base)
    n_float = log_epsilon / log_base
    n = math.ceil(n_float)

    # Print the equation with the numbers filled in, as requested.
    print("The formula to calculate n is: n = ceil(log(epsilon) / log(|1 - 2*delta|))")
    print("\nSubstituting the values:")
    print(f"n = ceil(log({epsilon}) / log(|1 - 2*{delta}|))")
    print(f"n = ceil({log_epsilon:.4f} / log(|{1 - 2 * delta:.2f}|))")
    print(f"n = ceil({log_epsilon:.4f} / log({base:.2f}))")
    print(f"n = ceil({log_epsilon:.4f} / {log_base:.4f})")
    print(f"n = ceil({n_float:.4f})")
    print(f"\nThe smallest integer n that satisfies the condition is: {n}")

solve_mixing_problem()
import math

def solve_mixing_problem(epsilon, delta):
    """
    Calculates the number of steps until the liquid in two cups is "the same".

    Args:
        epsilon (float): The maximum allowed difference in composition (0 < epsilon < 1).
        delta (float): The fraction of liquid exchanged in each step (0 < delta < 1).
    """
    print("--- Calculation ---")
    print(f"Goal: Find the smallest integer n such that the compositions differ by at most epsilon = {epsilon}.")
    print(f"Procedure: At each step, a delta = {delta} fraction of each cup is swapped.")
    print("")

    if not (0 < epsilon < 1 and 0 < delta < 1):
        print("Error: Epsilon and delta must be between 0 and 1.")
        return

    # Handle the special case where mixing is perfect in one step.
    if delta == 0.5:
        n = 1
        print("For delta = 0.5, perfect mixing occurs in a single step.")
        print(f"Final Answer: n = {n}")
        return

    # General case
    print("The derived formula for the number of steps n is:")
    print("n = ceil(log(epsilon) / log(abs(1 - 2*delta)))")
    print("\nPlugging in the values:")
    print(f"n = ceil(log({epsilon}) / log(abs(1 - 2*{delta})))")

    # Perform the calculations step by step
    numerator = math.log(epsilon)
    print(f"  Numerator: log({epsilon}) = {numerator:.4f}")

    base = abs(1 - 2 * delta)
    print(f"  Denominator base: abs(1 - 2*{delta}) = {base:.4f}")
    
    denominator = math.log(base)
    print(f"  Denominator: log({base:.4f}) = {denominator:.4f}")

    if denominator == 0:
        # This case is avoided by the delta constraints but good for robustness
        print("\nThe liquids will never mix sufficiently. The number of steps is infinite.")
        return
        
    ratio = numerator / denominator
    print(f"\nRatio: {numerator:.4f} / {denominator:.4f} = {ratio:.4f}")

    n = math.ceil(ratio)
    print(f"n = ceil({ratio:.4f})")
    
    print("\n--- Result ---")
    print(f"The procedure needs to be repeated {n} times.")


# --- Parameters ---
# You can change these values to see the result for different scenarios.
epsilon_val = 0.01  # Tolerance for "sameness"
delta_val = 0.1   # Fraction of liquid moved

# Execute the calculation
solve_mixing_problem(epsilon_val, delta_val)
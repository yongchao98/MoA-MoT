import math

def calculate_repetitions_and_explain(epsilon, delta):
    """
    Calculates the number of repetitions needed for the liquids to be "the same"
    and prints the derivation.

    Args:
        epsilon (float): The maximum allowed difference in concentration, between 0 and 1.
        delta (float): The fraction of liquid moved in each step, between 0 and 1.
    """
    if not (0 < epsilon < 1 and 0 < delta < 1):
        print("Error: Epsilon and delta must be between 0 and 1.")
        return

    print("The goal is to find the smallest integer n that satisfies the condition:")
    print("|concentration_A(n) - concentration_B(n)| <= epsilon")
    print("\nThis inequality can be simplified to: |1 - 2*delta|^n <= epsilon")

    # Handle the special case where mixing is perfect in one step.
    if delta == 0.5:
        print(f"\nFor delta = {delta}, mixing is perfect in a single step.")
        n = 1
        print("The number of repetitions needed is 1.")
        return

    # General case calculation
    print("\nTo solve for n, we take the natural logarithm of both sides:")
    print("n * ln(|1 - 2*delta|) <= ln(epsilon)")
    
    print("\nSince 0 < delta < 1 and delta != 0.5, |1 - 2*delta| is between 0 and 1.")
    print("This means ln(|1 - 2*delta|) is a negative number.")
    print("Therefore, when we divide by it, we must flip the inequality sign:")
    print("n >= ln(epsilon) / ln(|1 - 2*delta|)")

    abs_term = abs(1 - 2 * delta)
    log_epsilon = math.log(epsilon)
    log_abs_term = math.log(abs_term)
    
    n_float = log_epsilon / log_abs_term
    n_integer = math.ceil(n_float)

    print(f"\nSubstituting the given values epsilon = {epsilon} and delta = {delta}:")
    print(f"n >= ln({epsilon}) / ln(|1 - 2*{delta}|)")
    print(f"n >= ln({epsilon}) / ln({abs_term})")
    print(f"n >= {log_epsilon} / {log_abs_term}")
    print(f"n >= {n_float}")
    print(f"\nSince n must be an integer, we take the ceiling of the result.")
    print(f"The smallest integer n is {n_integer}.")

# --- Example Usage ---
# Define the parameters epsilon and delta.
# You can change these values to see the result for different parameters.
epsilon_val = 0.01
delta_val = 0.1

calculate_repetitions_and_explain(epsilon_val, delta_val)
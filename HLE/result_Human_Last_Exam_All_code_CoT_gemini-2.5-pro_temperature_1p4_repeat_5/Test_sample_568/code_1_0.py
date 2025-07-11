import math

def calculate_mixing_repetitions(epsilon, delta):
    """
    Calculates the number of repetitions needed for the liquid compositions
    in two cups to differ by at most epsilon.

    Args:
        epsilon (float): The maximum allowed difference, in (0, 1).
        delta (float): The fraction of liquid exchanged in each step, in (0, 1).
    """
    
    print(f"Calculating for epsilon = {epsilon} and delta = {delta}\n")
    
    # Handle edge cases where the general formula has a denominator of zero or is undefined.
    if delta <= 0 or delta >= 1:
        if epsilon < 1:
            print("Result: For delta <= 0 or delta >= 1, the liquids do not mix towards equilibrium.")
            print("The number of repetitions is effectively infinite.")
        else:
            print("Result: The condition is met from the start (n=0) as epsilon >= 1.")
        return

    if delta == 0.5:
        print("Special case: For delta = 0.5, perfect mixing is achieved in a single step.")
        print("Final Answer: n = 1")
        return

    # The general formula for the number of steps n is:
    # n = ceil( ln(epsilon) / ln(|1 - 2*delta|) )
    
    # Numerator of the argument to ceil
    log_epsilon = math.log(epsilon)
    
    # The base of the power in the inequality
    base = abs(1 - 2 * delta)
    
    # Denominator of the argument to ceil
    log_base = math.log(base)
    
    # The number of steps as a floating-point number
    n_float = log_epsilon / log_base
    
    # The final result must be an integer
    n_integer = math.ceil(n_float)

    print("The number of repetitions 'n' is the smallest integer satisfying:")
    print("n >= ln(epsilon) / ln(|1 - 2*delta|)\n")

    print("Step-by-step calculation:")
    print(f"n >= ln({epsilon}) / ln(|1 - 2*{delta}|)")
    print(f"n >= {log_epsilon:.4f} / ln(|{1 - 2 * delta}|)")
    print(f"n >= {log_epsilon:.4f} / ln({base:.4f})")
    print(f"n >= {log_epsilon:.4f} / {log_base:.4f}")
    print(f"n >= {n_float:.4f}\n")
    
    print(f"The smallest integer n that satisfies this is the ceiling of {n_float:.4f}.")
    print(f"Final Answer: n = {n_integer}")


# You can change these parameters to test other values.
# Both epsilon and delta must be between 0 and 1.
epsilon_param = 0.01
delta_param = 0.25

calculate_mixing_repetitions(epsilon_param, delta_param)

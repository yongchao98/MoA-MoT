import math

def solve_mixing_problem(epsilon, delta):
    """
    Calculates the number of times a procedure must be repeated for the
    liquid in two cups to be considered "the same".

    The formula derived is n = ceil(log(epsilon) / log(|1 - 2*delta|)).

    Args:
        epsilon (float): The tolerance for "sameness", must be in (0, 1).
        delta (float): The fraction of liquid exchanged, must be in (0, 1).
    """

    print(f"Starting calculation with parameters:\n  epsilon = {epsilon}\n  delta = {delta}")
    print("-" * 40)

    # Validate inputs based on problem constraints
    if not (0 < epsilon < 1 and 0 < delta < 1):
        print("Error: Epsilon and delta must be strictly between 0 and 1.")
        return

    # Handle the special case where delta = 0.5
    if delta == 0.5:
        n = 1
        print("This is a special case where delta = 0.5.")
        print("The liquids mix perfectly in a single step.")
        print(f"\nFinal Answer: The procedure must be repeated {n} time.")
        return

    # For the general case, we solve the inequality:
    # n >= log(epsilon) / log(|1 - 2*delta|)
    print("The required number of steps 'n' is the smallest integer satisfying:")
    print("n >= log(epsilon) / log(|1 - 2*delta|)")
    print("\nSubstituting the given values:")

    log_epsilon = math.log(epsilon)
    abs_term = abs(1 - 2 * delta)
    log_abs_term = math.log(abs_term)

    # Print the equation with numbers substituted
    print(f"n >= log({epsilon}) / log(|1 - 2*{delta}|)")
    print(f"n >= log({epsilon}) / log({abs_term:.4f})")
    print(f"n >= {log_epsilon:.4f} / {log_abs_term:.4f}")

    # Calculate the result
    n_float = log_epsilon / log_abs_term
    print(f"n >= {n_float:.4f}")

    # The number of steps must be an integer, so we take the ceiling.
    n = math.ceil(n_float)
    
    print(f"\nThe smallest integer 'n' that satisfies this is {n}.")
    print("-" * 40)
    print(f"Final Answer: The procedure must be repeated {n} times.")


# --- Example Usage ---
# You can change these values to see the result for different parameters.
epsilon_value = 0.01
delta_value = 0.1

solve_mixing_problem(epsilon_value, delta_value)
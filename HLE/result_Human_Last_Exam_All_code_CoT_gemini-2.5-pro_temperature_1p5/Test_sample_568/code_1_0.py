import math

def calculate_mixing_steps(epsilon, delta):
    """
    Calculates the number of steps to mix two liquids until their
    concentrations differ by at most epsilon.

    Args:
        epsilon (float): The maximum allowed difference in concentration, between 0 and 1.
        delta (float): The fraction of liquid exchanged in each step, between 0 and 1.
    """
    # You can change the input parameters here
    # epsilon = 0.01
    # delta = 0.1

    if not (0 < epsilon < 1 and 0 < delta < 1):
        print("Error: Epsilon and delta must be in the interval (0, 1).")
        return

    print(f"Solving for epsilon = {epsilon} and delta = {delta}")
    print("-" * 30)

    # Handle the special case where delta = 0.5
    if delta == 0.5:
        n = 1
        print("When delta = 0.5, the liquids mix perfectly in one step.")
        print("The difference in concentrations becomes 0, which is <= epsilon.")
        print(f"The number of repetitions is: {n}")
        return

    # For delta != 0.5, we solve the inequality: |(1 - 2*delta)^n| <= epsilon
    # By taking logs: n * log(|1 - 2*delta|) <= log(epsilon)
    # Then: n >= log(epsilon) / log(|1 - 2*delta|)

    term = 1 - 2 * delta
    abs_term = math.fabs(term)
    log_epsilon = math.log(epsilon)
    log_abs_term = math.log(abs_term)

    print("The number of steps n must satisfy the inequality:")
    print(f"|1 - 2 * ({delta})|^n <= {epsilon}")
    print(f"|{term:.2f}|^n <= {epsilon}")
    print(f"{abs_term:.4f}^n <= {epsilon}")
    print("\nTaking the natural logarithm of both sides:")
    print(f"n * log({abs_term:.4f}) <= log({epsilon})")
    print(f"n * ({log_abs_term:.4f}) <= {log_epsilon:.4f}")

    print("\nBecause log(|1 - 2*delta|) is negative, we reverse the inequality sign when dividing:")
    n_float = log_epsilon / log_abs_term
    print(f"n >= {log_epsilon:.4f} / {log_abs_term:.4f}")
    print(f"n >= {n_float:.4f}")

    # The number of repetitions must be an integer, so we take the ceiling.
    n = math.ceil(n_float)

    print("\nSince n must be an integer, we take the ceiling of the result:")
    print(f"n = ceil({n_float:.4f})")
    print(f"\nThe final number of repetitions required is: {int(n)}")

# --- User-defined parameters ---
# You can change these values to see the result for different scenarios.
epsilon_value = 0.01
delta_value = 0.25
# --- End of user-defined parameters ---

calculate_mixing_steps(epsilon_value, delta_value)
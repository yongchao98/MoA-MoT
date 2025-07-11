import math

def calculate_mixing_repetitions(epsilon, delta):
    """
    Calculates the number of repetitions (n) needed for the liquid in two cups
    to be considered "the same" based on the parameters epsilon and delta.

    Args:
        epsilon (float): The maximum allowed difference in content fraction (0 < epsilon < 1).
        delta (float): The fraction of liquid exchanged in each step (0 < delta < 1).
    """

    # --- Step-by-step derivation explanation ---
    print("--- Derivation of the formula ---")
    print("Let r_A(n) be the fraction of red liquid in cup A after n steps.")
    print("Let r_B(n) be the fraction of red liquid in cup B after n steps.")
    print("Initially, r_A(0) = 1 and r_B(0) = 0.")
    print("\nIn each step, we exchange a delta-fraction of liquid between the cups.")
    print("The amount of red liquid in cup A for the next step, r_A(n+1), is:")
    print("r_A(n+1) = r_A(n) - delta*r_A(n) + delta*r_B(n)")
    print("Since the total red liquid is 1, r_A(n) + r_B(n) = 1. So, r_B(n) = 1 - r_A(n).")
    print("r_A(n+1) = r_A(n) - delta*r_A(n) + delta*(1 - r_A(n))")
    print("r_A(n+1) = (1 - 2*delta)*r_A(n) + delta")
    print("\nSolving this recurrence relation gives the explicit formula for r_A(n):")
    print("r_A(n) = 0.5 + 0.5 * (1 - 2*delta)^n")
    print("\nThe difference in red liquid concentration between the cups is |r_A(n) - r_B(n)|.")
    print("|r_A(n) - r_B(n)| = |(0.5 + 0.5*(1-2d)^n) - (0.5 - 0.5*(1-2d)^n)| = |(1 - 2*delta)^n|")
    print("\nWe need to find the smallest integer n where this difference is at most epsilon:")
    print("|1 - 2*delta|^n <= epsilon")
    print("\nTaking the logarithm and solving for n, we get:")
    print("n >= log(epsilon) / log(|1 - 2*delta|)")
    print("Since n must be an integer, we take the ceiling.")
    print("n = ceil(log(epsilon) / log(|1 - 2*delta|))")
    print("---------------------------------------\n")

    # --- Calculation for the given parameters ---
    print(f"Given parameters:")
    print(f"epsilon = {epsilon}")
    print(f"delta = {delta}\n")

    if not (0 < epsilon < 1 and 0 < delta < 1):
        print("Error: Epsilon and delta must be between 0 and 1.")
        return

    # Handle the special case where delta = 0.5, which leads to division by zero in the formula.
    if delta == 0.5:
        # After 1 step, the cups are perfectly mixed. The difference is 0, which is <= epsilon.
        n = 1
        print("For delta = 0.5, the liquids mix perfectly in one step.")
        print("n = 1")
    else:
        # General case calculation
        log_epsilon = math.log(epsilon)
        log_abs_term = math.log(abs(1 - 2 * delta))
        result_float = log_epsilon / log_abs_term
        n = math.ceil(result_float)

        print("The final equation with the given numbers is:")
        print(f"n = ceil(log({epsilon}) / log(|1 - 2 * {delta}|))")
        print(f"n = ceil({log_epsilon} / {log_abs_term})")
        print(f"n = ceil({result_float})")
        print(f"n = {int(n)}")

    print("\nThus, the procedure must be repeated", int(n), "times.")


# You can change these parameters to see the result for different values.
# Fix parameters epsilon, delta in (0, 1)
epsilon_param = 0.01
delta_param = 0.1

calculate_mixing_repetitions(epsilon_param, delta_param)
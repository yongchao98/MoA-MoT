import math

def solve_mixing_problem(epsilon, delta):
    """
    Calculates the number of repetitions needed for the liquids in two cups to be "the same".

    Args:
        epsilon (float): The maximum allowed difference in content fractions, between 0 and 1.
        delta (float): The fraction of liquid exchanged in each step, between 0 and 1.
    """
    # Input validation
    if not (0 < epsilon < 1 and 0 < delta < 1):
        print("Error: epsilon and delta must be between 0 and 1.")
        return

    print("This script calculates the number of steps 't' required for the liquid concentrations to be within epsilon of each other.")
    print(f"Given parameters: epsilon = {epsilon}, delta = {delta}\n")

    # The recurrence relation for the fraction of red liquid in cup A, r_A(t), is:
    # r_A(t+1) = r_A(t) * (1 - 2*delta) + delta
    # The solution to this is r_A(t) = 0.5 + 0.5 * (1 - 2*delta)^t.
    #
    # The condition for the liquids to be "the same" is |r_A(t) - r_B(t)| <= epsilon,
    # which simplifies to |(1 - 2*delta)^t| <= epsilon.
    #
    # Solving for t gives: t >= log(epsilon) / log(|1 - 2*delta|)
    # Since t must be an integer, we take the ceiling.

    print("The required number of steps 't' is the smallest integer satisfying the inequality.")
    print("The derived formula is: t = ceil( log(epsilon) / log(|1 - 2*delta|) )\n")

    # Handle the special case where delta = 0.5
    if delta == 0.5:
        # If delta is 0.5, the cups are perfectly mixed in one step.
        # The denominator log(|1 - 2*0.5|) = log(0) would be mathematically undefined.
        t = 1
        print("Calculation:")
        print("For delta = 0.5, perfect mixing is achieved in a single step.")
        print(f"Therefore, t = {t}")
    else:
        # Perform the calculation for the general case
        abs_val = abs(1 - 2 * delta)
        log_epsilon = math.log(epsilon)
        log_abs_val = math.log(abs_val)
        
        # Check for delta being 0 or 1, which are excluded by problem constraints but good to handle
        if log_abs_val == 0:
            print("For delta = 0 or 1, the liquids never mix to the desired level.")
            print("The number of steps is infinite.")
            return

        result = log_epsilon / log_abs_val
        t = math.ceil(result)

        print("Calculation with the given numbers:")
        print(f"t = ceil( log({epsilon}) / log(|1 - 2 * {delta}|) )")
        print(f"t = ceil( {log_epsilon:.4f} / log(|{1 - 2 * delta}|) )")
        print(f"t = ceil( {log_epsilon:.4f} / {log_abs_val:.4f} )")
        print(f"t = ceil( {result:.4f} )")
        print(f"t = {t}")

    print(f"\nConclusion: The procedure must be repeated {t} times.")

# --- User-configurable parameters ---
# ε: The maximum difference in content fractions for the liquids to be considered "the same".
epsilon_param = 0.01

# δ: The fraction of liquid to be moved from each cup.
delta_param = 0.25
# ------------------------------------

solve_mixing_problem(epsilon_param, delta_param)
import math

def solve_mixing_problem():
    """
    This script calculates the number of steps for the liquids in two cups
    to be considered 'the same' based on the derived formula.
    """
    # Fix parameters for the demonstration, as per the problem description.
    epsilon = 0.01
    delta = 0.1

    print("This script calculates the number of steps 'n' for the liquid mixing problem.")
    print(f"We are given epsilon = {epsilon} and delta = {delta}.\n")

    print("The formula to find n is: n = ceil(log(epsilon) / log(|1 - 2*delta|))\n")

    # Handle the special case where delta = 0.5, which would cause a math error.
    if delta == 0.5:
        # In this case, perfect mixing is achieved in one step.
        n = 1
        print("Special case: For delta = 0.5, perfect mixing occurs in 1 step.")
    else:
        # Calculate the result using the general formula
        log_epsilon = math.log(epsilon)
        base = abs(1 - 2 * delta)
        log_base = math.log(base)
        raw_n = log_epsilon / log_base
        n = math.ceil(raw_n)

        print("Step 1: Calculate the terms of the formula.")
        print(f"  epsilon           = {epsilon}")
        print(f"  delta             = {delta}")
        print(f"  log(epsilon)      = log({epsilon}) = {log_epsilon:.4f}")
        print(f"  |1 - 2*delta|     = |1 - 2*{delta}| = {base:.4f}")
        print(f"  log(|1 - 2*delta|) = log({base:.4f}) = {log_base:.4f}\n")

        print("Step 2: Plug the values into the inequality n >= log(epsilon) / log(|1 - 2*delta|).")
        print(f"  n >= {log_epsilon:.4f} / {log_base:.4f}")
        print(f"  n >= {raw_n:.4f}\n")

        print(f"Step 3: The smallest integer n is the ceiling of {raw_n:.4f}.")
        print(f"  n = {n}\n")

    print("--- Result ---")
    print(f"The procedure must be repeated {n} times for the liquids to be considered 'the same'.")

solve_mixing_problem()
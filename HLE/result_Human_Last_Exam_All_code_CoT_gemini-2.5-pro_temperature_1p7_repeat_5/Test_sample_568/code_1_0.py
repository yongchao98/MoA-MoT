import math

def solve_mixing_problem():
    """
    Calculates the number of steps required for the liquids in two cups
    to be considered "the same" based on the mixing procedure described.
    """
    
    # Parameters from the problem.
    # You can change these values to see the result for different conditions.
    # epsilon is the maximum allowed difference in concentration (must be in (0, 1)).
    # delta is the fraction of liquid exchanged in each step (must be in (0, 1)).
    epsilon = 0.01
    delta = 0.1

    print(f"Calculating the number of steps 'n' for epsilon = {epsilon} and delta = {delta}")
    print("---------------------------------------")

    # Validate inputs
    if not (0 < epsilon < 1 and 0 < delta < 1):
        print("Error: Epsilon and delta must be between 0 and 1.")
        return

    # Handle the special case where delta = 0.5
    if delta == 0.5:
        # If delta is 0.5, the cups are perfectly mixed in one step.
        # The concentration of red liquid becomes 0.5 in both cups,
        # so the difference is 0, which is less than any positive epsilon.
        n = 1
        print("Special case: delta = 0.5")
        print("The liquids mix perfectly in one step.")

    else:
        # The general formula is derived from the inequality: |1 - 2*delta|^n <= epsilon
        # Solving for n gives: n >= log(epsilon) / log(abs(1 - 2*delta))
        print("The number of repetitions 'n' is the smallest integer satisfying the inequality:")
        print("  |1 - 2*delta|^n  <=  epsilon")
        print("\nSolving for n, we get:")
        print("  n  >=  log(epsilon) / log(abs(1 - 2*delta))")
        print("Since n must be an integer, we take the ceiling of the result.")
        
        print("\nSubstituting the given values:")
        
        val_in_abs = 1 - 2 * delta
        abs_val = abs(val_in_abs)
        log_epsilon = math.log(epsilon)
        log_abs_val = math.log(abs_val)
        
        print(f"  n  >=  log({epsilon}) / log(abs(1 - 2*{delta}))")
        print(f"  n  >=  log({epsilon}) / log(abs({val_in_abs:.2f}))")
        print(f"  n  >=  log({epsilon}) / log({abs_val})")
        print(f"  n  >=  {log_epsilon:.6f} / {log_abs_val:.6f}")
        
        ratio = log_epsilon / log_abs_val
        print(f"  n  >=  {ratio:.6f}")
        
        n = math.ceil(ratio)

    print("\n=======================================")
    print(f"The final number of repetitions required is: {n}")
    print("=======================================")

# Execute the function
solve_mixing_problem()
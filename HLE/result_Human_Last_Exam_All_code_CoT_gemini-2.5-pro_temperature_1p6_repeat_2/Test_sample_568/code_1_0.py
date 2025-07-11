import math

def solve_mixing_problem(epsilon, delta):
    """
    Calculates the number of steps until the liquids in two cups are "the same".

    Args:
        epsilon (float): The tolerance for the difference in liquid composition, must be in (0, 1).
        delta (float): The fraction of liquid exchanged in each step, must be in (0, 1).
    """
    
    if not (0 < epsilon < 1):
        print("Error: Epsilon (epsilon) must be between 0 and 1.")
        return
    if not (0 < delta < 1):
        print("Error: Delta (delta) must be between 0 and 1.")
        return

    print(f"Given parameters:\nepsilon = {epsilon}\ndelta = {delta}\n")
    print("The minimum number of steps 'n' is the smallest integer satisfying the equation:")
    print("|(1 - 2*delta)^n| <= epsilon\n")
    print("This can be solved for n as: n >= log(epsilon) / log(|1 - 2*delta|)")
    print("So, n = ceil(log(epsilon) / log(|1 - 2*delta|))\n")
    
    # Handle the special case where delta = 0.5
    if delta == 0.5:
        # After 1 step, the mixture is perfect (50/50), so the difference is 0.
        # Since epsilon > 0, the condition 0 <= epsilon is met.
        n = 1
        print("For the special case where delta = 0.5, the cups are perfectly mixed in 1 step.")
        final_result = 1
    else:
        # General case
        print("Calculation steps:")
        
        # Denominator of the formula
        log_arg = abs(1 - 2 * delta)
        
        # Numerator value
        log_epsilon = math.log(epsilon)
        
        # Denominator value
        log_abs_val = math.log(log_arg)
        
        n_float = log_epsilon / log_abs_val
        
        # Ceiling to get the final integer number of steps
        n_int = math.ceil(n_float)

        print(f"n = ceil(log({epsilon}) / log(|1 - 2*{delta}|))")
        print(f"n = ceil({log_epsilon:.4f} / log(|{1 - 2 * delta:.4f}|))")
        print(f"n = ceil({log_epsilon:.4f} / log({log_arg:.4f}))")
        print(f"n = ceil({log_epsilon:.4f} / {log_abs_val:.4f})")
        print(f"n = ceil({n_float:.4f})")
        final_result = n_int

    print("\n----------------------------------")
    print(f"Final Answer: The procedure needs to be repeated {final_result} times.")
    print("----------------------------------")
    
    # Returning the final answer as per the problem format requirement
    return final_result


# Example usage with some fixed parameters
# You can change these values to test other scenarios
epsilon_val = 0.01
delta_val = 0.1

final_answer = solve_mixing_problem(epsilon_val, delta_val)
# The format <<<answer>>> is for programmatic evaluation
# and contains the final numerical result.
print(f"\n<<<{final_answer}>>>")

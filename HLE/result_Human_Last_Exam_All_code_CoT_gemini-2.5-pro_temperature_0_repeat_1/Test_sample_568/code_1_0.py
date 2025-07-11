import math

def solve_mixing_problem(epsilon, delta):
    """
    Calculates the number of steps required for the liquid in two cups to become "the same".

    Args:
        epsilon (float): The maximum allowed difference in composition (0 < epsilon < 1).
        delta (float): The fraction of liquid exchanged in each step (0 < delta < 1).
    """
    print(f"The user has provided the following parameters: epsilon = {epsilon}, delta = {delta}\n")

    print("The problem is to find the number of repetitions, t, until the liquids are the same.")
    print("The condition for the liquids being 'the same' is when the difference in the concentration of red liquid in cup A and cup B is at most epsilon.")
    print("This can be expressed with the inequality: |1 - 2*delta|^t <= epsilon\n")

    print("To solve for t, we take the natural logarithm of both sides:")
    print("t * ln(|1 - 2*delta|) <= ln(epsilon)")
    print("Since |1 - 2*delta| < 1, its logarithm is negative, so we reverse the inequality when dividing:")
    print("t >= ln(epsilon) / ln(|1 - 2*delta|)\n")
    
    print("Since t must be an integer, we take the ceiling of the result.")
    print("The final formula is: t = ceil(ln(epsilon) / ln(|1 - 2*delta|))\n")

    # Handle the special case where delta = 0.5
    if delta == 0.5:
        t = 1
        print("For the special case where delta = 0.5, perfect mixing occurs in one step.")
        print(f"The number of repetitions t is: {t}")
        return

    # Perform the calculation
    log_epsilon = math.log(epsilon)
    abs_val = abs(1 - 2 * delta)
    log_abs_val = math.log(abs_val)
    
    # Check for invalid inputs that might lead to math errors, though problem constraints should prevent this.
    if log_abs_val == 0:
        print("Error: delta cannot be 0 or 1.")
        return
    if epsilon <= 0 or epsilon >= 1:
        print("Error: epsilon must be between 0 and 1.")
        return

    result_float = log_epsilon / log_abs_val
    t = math.ceil(result_float)

    print("Now, let's plug in the numbers into the formula:")
    # Final equation with numbers
    print(f"t = ceil(ln({epsilon}) / ln(|1 - 2*{delta}|))")
    print(f"t = ceil({log_epsilon} / ln(|{1 - 2*delta}|))")
    print(f"t = ceil({log_epsilon} / {log_abs_val})")
    print(f"t = ceil({result_float})")
    print(f"The final result for the number of repetitions t is: {t}")


# --- Parameters ---
# You can change these values to solve for different epsilon and delta.
epsilon = 0.01  # The "sameness" threshold
delta = 0.1     # The fraction of liquid to exchange

# Run the solver
solve_mixing_problem(epsilon, delta)

# The final answer is the formula derived and implemented.
# For the given parameters epsilon=0.01 and delta=0.1, the answer is 21.
# The general answer is the function itself.
# Let's express the final answer as the formula string.
final_answer_formula = "ceil(ln(epsilon) / ln(|1 - 2*delta|))"
# For the specific values, the answer is an integer.
final_answer_value = math.ceil(math.log(epsilon) / math.log(abs(1 - 2 * delta))) if delta != 0.5 else 1

# The prompt asks for a single answer at the end. The most general answer is the formula.
# However, since the code calculates a specific value, let's output that value.
# For the default parameters epsilon=0.01, delta=0.1, the result is 21.
print(f"\n<<< {final_answer_value} >>>")
import math

def solve_mixing_problem(epsilon, delta):
    """
    Calculates the number of repetitions (n) needed for the liquid in two cups
    to be considered "the same" based on the given parameters epsilon and delta.

    The state is tracked by r_a(n), the fraction of red liquid in cup A after n steps.
    Initially, r_a(0) = 1. The recurrence relation derived from the mixing
    procedure is r_a(n+1) = (1 - 2*delta) * r_a(n) + delta.
    The solution to this recurrence is r_a(n) = 0.5 + 0.5 * (1 - 2*delta)^n.

    The "sameness" condition is |r_a(n) - r_b(n)| <= epsilon, where r_b(n) = 1 - r_a(n).
    This simplifies to |2*r_a(n) - 1| <= epsilon.
    Substituting the solution for r_a(n), we get: |(1 - 2*delta)^n| <= epsilon,
    which is equivalent to |1 - 2*delta|^n <= epsilon.

    Solving for n gives: n >= log(epsilon) / log(|1 - 2*delta|).
    Since n must be an integer, we take the ceiling of the result.
    """

    print(f"Given parameters: epsilon = {epsilon}, delta = {delta}\n")
    print("The formula to calculate the number of repetitions (n) is:")
    print("n = ceil(log(epsilon) / log(|1 - 2*delta|))\n")

    if not (0 < epsilon < 1 and 0 < delta < 1):
        print("Error: Epsilon and Delta must be in the range (0, 1).")
        return
        
    # Handle the special case where delta = 0.5
    if delta == 0.5:
        n = 1
        print("For delta = 0.5, the denominator log(|1 - 2*0.5|) = log(0) is undefined.")
        print("In this case, the liquids mix perfectly in a single step.")
        print(f"The number of repetitions needed is n = {n}.")
        # The following line is for the final answer block
        global final_answer
        final_answer = n
        return

    # Calculate n for the general case
    abs_val = abs(1 - 2 * delta)
    log_epsilon = math.log(epsilon)
    log_abs_val = math.log(abs_val)
    
    # Ensure we don't divide by zero if abs_val is 1 (delta=0 or delta=1)
    # Though the problem statement implies 0 < delta < 1
    if log_abs_val == 0:
        print("For this value of delta, the liquids never mix to the desired degree.")
        print("The number of repetitions would be infinite.")
        return

    ratio = log_epsilon / log_abs_val
    n = math.ceil(ratio)
    
    # The following line is for the final answer block
    final_answer = n

    print("Plugging in the values:")
    print(f"n = ceil(log({epsilon}) / log(|1 - 2*{delta}|))")
    print(f"n = ceil(log({epsilon}) / log({abs_val:.4f}))")
    print(f"n = ceil({log_epsilon:.4f} / {log_abs_val:.4f})")
    print(f"n = ceil({ratio:.4f})")
    print(f"n = {n}")


# Example parameters, as none were specified in the prompt.
# You can change these values to solve for different parameters.
epsilon_val = 0.01
delta_val = 0.1

# This variable will be captured for the final answer
final_answer = None
solve_mixing_problem(epsilon_val, delta_val)
print(f"\n<<<>>>\n{final_answer}")
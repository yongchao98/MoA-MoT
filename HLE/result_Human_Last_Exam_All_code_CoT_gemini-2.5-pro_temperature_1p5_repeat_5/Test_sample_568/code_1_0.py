import math

def solve_mixing_problem(epsilon, delta):
    """
    Calculates the number of repetitions (n) needed for the two cups
    to be considered "the same" based on the given epsilon and delta.

    The state of the system can be described by the fraction of red liquid
    in cup A, R_A(n). This follows the recurrence:
    R_A(n+1) = (1 - 2*delta) * R_A(n) + delta
    
    The closed-form solution is:
    R_A(n) = 0.5 + 0.5 * (1 - 2*delta)^n
    
    The difference in red liquid fraction between the cups is |R_A(n) - R_B(n)|,
    which simplifies to |(1 - 2*delta)^n|.
    
    The condition to be "the same" is |(1 - 2*delta)^n| <= epsilon.
    
    Solving for n gives: n >= log(epsilon) / log(|1 - 2*delta|).
    Since n must be an integer, we take the ceiling of this value.
    """

    print(f"Calculating for epsilon = {epsilon} and delta = {delta}\n")

    # The parameters epsilon and delta must be in the interval (0, 1).
    if not (0 < epsilon < 1 and 0 < delta < 1):
        print("Error: Epsilon and delta must be between 0 and 1.")
        return

    # Handle the special case where mixing is perfect in one step.
    if delta == 0.5:
        n = 1
        print("When delta is 0.5, the liquids mix perfectly in one step.")
        print(f"Final answer: n = {n}")
        return n

    # The equation to solve for n:
    # n >= log(epsilon) / log(|1 - 2*delta|)
    
    log_epsilon = math.log(epsilon)
    base = abs(1 - 2 * delta)
    log_base = math.log(base)
    
    ratio = log_epsilon / log_base
    n = math.ceil(ratio)

    print("The required number of repetitions 'n' is the smallest integer satisfying:")
    print(f"|1 - 2*\u03B4|^n \u2264 \u03B5")
    print(f"|1 - 2*{delta}|^{n} \u2264 {epsilon}\n")
    
    print("Taking the logarithm of both sides and solving for n:")
    print(f"n \u2265 log(\u03B5) / log(|1 - 2*\u03B4|)")
    print(f"n \u2265 log({epsilon}) / log(|1 - 2*{delta}|)")
    print(f"n \u2265 {log_epsilon:.4f} / log({base:.4f})")
    print(f"n \u2265 {log_epsilon:.4f} / {log_base:.4f}")
    print(f"n \u2265 {ratio:.4f}\n")
    
    print(f"The smallest integer n is the ceiling of {ratio:.4f}.")
    print(f"Final answer: n = {n}")
    
    return n

# --- Main execution ---
# Fix parameters epsilon and delta as per the problem description.
# Let's use epsilon = 0.01 (1% difference) and delta = 0.25 (25% fraction).
epsilon_val = 0.01
delta_val = 0.25

solve_mixing_problem(epsilon_val, delta_val)
import math

def solve_mixing_problem(epsilon, delta):
    """
    Calculates the number of steps until the liquid in two cups is "the same".

    Args:
      epsilon (float): The tolerance for the difference in liquid composition, between 0 and 1.
      delta (float): The fraction of liquid moved in each step, between 0 and 1.
    """
    if not (0 < epsilon < 1 and 0 < delta < 1):
        print("Error: Epsilon and delta must be between 0 and 1.")
        return

    print(f"Given parameters: epsilon = {epsilon}, delta = {delta}")
    print("\nStep 1: The recurrence relation for the fraction of red liquid in cup A, r_A(n), is:")
    print("r_A(n) = r_A(n-1) * (1 - 2*delta) + delta")
    print("\nStep 2: The solved closed-form for r_A(n) is:")
    print("r_A(n) = 0.5 + 0.5 * (1 - 2*delta)^n")
    print("\nStep 3: The condition for the liquids to be 'the same' is |r_A(n) - r_B(n)| <= epsilon, which simplifies to:")
    print("|(1 - 2*delta)^n| <= epsilon")
    print("\nStep 4: Solving for n by taking logarithms gives:")
    print("n >= log(epsilon) / log(|1 - 2*delta|)")
    print("\nSince n must be an integer, we take the ceiling of the result.")

    # Handle the special case where delta = 0.5
    if delta == 0.5:
        n = 1
        print("\nCalculation:")
        print("When delta = 0.5, the liquids mix perfectly in one step.")
        print(f"Final Answer (n) = {n}")
        return n

    # Calculation for the general case
    abs_term = abs(1 - 2 * delta)
    log_epsilon = math.log(epsilon)
    log_delta_term = math.log(abs_term)
    
    # Avoid division by zero if abs_term is 1 (i.e., delta is 0 or 1, which is ruled out)
    if log_delta_term == 0:
        print("\nWith the given delta, the mixture never changes. n is infinite.")
        return float('inf')

    result_frac = log_epsilon / log_delta_term
    n = math.ceil(result_frac)

    print("\nCalculation with the given numbers:")
    final_equation = f"n = ceil(log({epsilon}) / log(|1 - 2*{delta}|))"
    print(final_equation)
    final_equation_eval = f"n = ceil({log_epsilon:.4f} / log({abs_term:.4f}))"
    print(final_equation_eval)
    final_equation_res = f"n = ceil({result_frac:.4f})"
    print(final_equation_res)
    print(f"Final Answer (n) = {n}")
    
    return n

# --- Example Usage ---
# You can change these values to see the result for different parameters.
epsilon_val = 0.01
delta_val = 0.1

# Run the solver and get the final numerical answer
final_n = solve_mixing_problem(epsilon_val, delta_val)
print(f"\n<<< {final_n} >>>")
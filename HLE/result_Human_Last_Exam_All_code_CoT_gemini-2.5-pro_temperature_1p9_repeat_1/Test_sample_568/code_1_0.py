import math

def solve_mixing_problem(epsilon, delta):
    """
    Calculates the number of repetitions needed for the liquids in two cups
    to be considered "the same".

    Args:
        epsilon (float): The tolerance for the difference, between 0 and 1.
        delta (float): The fraction of liquid moved, between 0 and 1.
    """
    print(f"Given parameters ε = {epsilon} and δ = {delta}\n")
    
    if not (0 < epsilon < 1 and 0 < delta < 1):
        print("Error: Epsilon (ε) and delta (δ) must be in the interval (0, 1).")
        return

    # The final answer for n is given by the formula:
    # n = 1, if δ = 1/2
    # n = ceil(ln(ε) / ln(|1 - 2δ|)), if δ ≠ 1/2
    # We will print the equation and the values used for the calculation.

    print("The number of repetitions 'n' must satisfy the inequality:")
    print("|1 - 2*δ|^n <= ε")
    print("\nSolving for n, we get:")
    print("n >= ln(ε) / ln(|1 - 2*δ|)\n")

    if delta == 0.5:
        n = 1
        print("Since δ = 0.5, the term |1 - 2*δ| is 0.")
        print("The system reaches equilibrium in a single step.")
        print(f"Final Answer (n): {n}")
        return n

    # Calculate the terms in the equation
    log_epsilon = math.log(epsilon)
    abs_term = abs(1 - 2 * delta)
    log_abs_term = math.log(abs_term)
    
    # Check for invalid delta that makes ln(|1-2d|) zero.
    # This only happens if |1-2d|=1, which means d=0 or d=1, but we already constrained d.
    if log_abs_term == 0:
        print("Result is undefined as the mixture never changes (δ=0 or δ=1).")
        return

    ratio = log_epsilon / log_abs_term
    n = math.ceil(ratio)

    print("Substituting the given values into the equation:")
    # Using the variable names for clarity in the equation printout
    final_equation = f"n >= ln({epsilon}) / ln(|1 - 2*{delta}|)"
    print(final_equation)
    
    # Printing each number in the calculation
    calc_step1 = f"n >= {log_epsilon} / ln({abs_term})"
    print(calc_step1)

    calc_step2 = f"n >= {log_epsilon} / {log_abs_term}"
    print(calc_step2)

    calc_step3 = f"n >= {ratio}"
    print(calc_step3)
    
    print("\nSince n must be an integer, we take the ceiling of the result.")
    print(f"Final Answer (n): {n}")
    return n

# Example usage with sample parameters, as they were not provided in the prompt.
# You can change these values to see the result for different parameters.
epsilon_value = 0.01
delta_value = 0.1

solve_mixing_problem(epsilon_value, delta_value)
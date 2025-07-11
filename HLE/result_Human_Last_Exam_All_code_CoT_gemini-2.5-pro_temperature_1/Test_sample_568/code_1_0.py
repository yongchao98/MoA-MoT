import math

def solve_mixing_problem(epsilon, delta):
    """
    Calculates the number of repetitions needed for the liquids to be 'the same'.

    Args:
        epsilon (float): The tolerance for the difference in composition, must be in (0, 1).
        delta (float): The fraction of liquid moved in each step, must be in (0, 1).

    Returns:
        int: The minimum number of repetitions N.
    """
    # Validate input parameters based on the problem description
    if not (0 < epsilon < 1 and 0 < delta < 1):
        raise ValueError("Epsilon and delta must be in the range (0, 1).")

    # Handle the special case where delta = 0.5
    # Equilibrium is reached in a single step.
    if delta == 0.5:
        return 1

    # For other cases, calculate N using the derived formula.
    # We need to find the smallest integer N such that: |1 - 2*delta|^N <= epsilon
    # This leads to N >= log(epsilon) / log(|1 - 2*delta|)
    
    # The base of the logarithm does not matter as long as it's consistent.
    # Python's math.log is the natural logarithm (base e).
    log_epsilon = math.log(epsilon)
    log_abs_val = math.log(abs(1 - 2 * delta))
    
    # The result is the ceiling of the division.
    N_float = log_epsilon / log_abs_val
    N = math.ceil(N_float)
    
    return int(N)

# --- Main execution block ---
# You can change these parameters to test different scenarios.
epsilon_param = 0.01
delta_param = 0.1

try:
    # Calculate the result
    repetitions = solve_mixing_problem(epsilon_param, delta_param)

    # Print the output
    print(f"Given parameters:")
    print(f"epsilon = {epsilon_param}")
    print(f"delta = {delta_param}")
    print("-" * 30)
    print(f"The minimum number of repetitions required is: {repetitions}")
    print("-" * 30)
    print("This result is the smallest integer N satisfying the equation:")
    
    # Fulfilling the request to "output each number in the final equation"
    if delta_param == 0.5:
        print("N = 1 (since delta = 0.5, equilibrium is reached in one step)")
    else:
        print(f"N = ceil(log(epsilon) / log(|1 - 2*delta|))")
        print(f"N = ceil(log({epsilon_param}) / log(|1 - 2*{delta_param}|))")
        log_eps = math.log(epsilon_param)
        log_abs = math.log(abs(1 - 2*delta_param))
        print(f"N = ceil({log_eps:.4f} / {log_abs:.4f})")
        print(f"N = ceil({log_eps / log_abs:.4f})")
        print(f"N = {repetitions}")

except ValueError as e:
    print(f"An error occurred: {e}")

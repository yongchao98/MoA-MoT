import math

def calculate_mixing_repetitions(epsilon, delta):
    """
    Calculates the number of repetitions (n) needed for the liquid concentrations
    in two cups to differ by at most epsilon.

    Args:
        epsilon (float): The maximum allowed difference in concentration, between 0 and 1.
        delta (float): The fraction of liquid exchanged in each step, between 0 and 1.
    """
    print(f"Parameters: ε = {epsilon}, δ = {delta}\n")

    # Handle edge case where the condition is already met.
    # The initial difference is 1. If epsilon is 1 or greater, 0 steps are needed.
    if epsilon >= 1:
        n = 0
        print("The initial difference is 1. Since ε >= 1, the condition is met from the start.")
        print(f"\nNumber of repetitions needed: {n}")
        return

    # Handle edge case where mixing is perfect in one step.
    if delta == 0.5:
        n = 1
        print("When δ = 0.5, the liquids mix perfectly in a single step.")
        print(f"The difference becomes 0, which is less than or equal to any positive ε.")
        print(f"\nNumber of repetitions needed: {n}")
        return

    # Calculate the base of the exponent, which is |1 - 2δ|.
    base = abs(1 - 2 * delta)

    # Explain the equation to be solved: |1 - 2δ|^n ≤ ε
    print("We need to find the smallest integer n that satisfies the inequality:")
    print(f"|1 - 2*δ|^n <= ε")
    print("Using the provided values, the equation is:")
    print(f"{base:.4f}^n <= {epsilon}")
    
    # Solve for n by taking logarithms: n ≥ log(ε) / log(|1 - 2δ|)
    print("\nTo solve for n, we use logarithms:")
    print(f"n >= log({epsilon}) / log({base:.4f})")
    
    # Perform the calculation
    log_epsilon = math.log(epsilon)
    log_base = math.log(base)
    result = log_epsilon / log_base
    
    print(f"n >= {log_epsilon:.4f} / {log_base:.4f}")
    print(f"n >= {result:.4f}")
    
    # The number of steps must be an integer, so we take the ceiling.
    n = math.ceil(result)
    
    print(f"\nThe smallest integer n that satisfies this is the ceiling of {result:.4f}.")
    print(f"\nFinal Answer: The procedure must be repeated {n} times.")


# --- Main execution ---
# Set the parameters for the problem.
# You can change these values to see different outcomes.
epsilon_param = 0.01  # We want the concentrations to differ by at most 1%.
delta_param = 0.1   # In each step, 10% of the liquid is exchanged.

calculate_mixing_repetitions(epsilon_param, delta_param)
import math

def solve_mixing_problem(epsilon, delta):
    """
    Calculates the number of repetitions 'n' for the liquid in two cups to be considered "the same".

    Args:
        epsilon (float): The maximum allowed difference in concentration (must be in (0, 1)).
        delta (float): The fraction of liquid exchanged in each step (must be in (0, 1)).
    """
    print(f"This script calculates the number of repetitions 'n' required.")
    print(f"Parameters: ε = {epsilon}, δ = {delta}\n")

    # Validate inputs
    if not (0 < epsilon < 1 and 0 < delta < 1):
        print("Error: Epsilon (ε) and delta (δ) must both be in the interval (0, 1).")
        return

    # Handle the special case where delta = 0.5
    if delta == 0.5:
        n = 1
        print("For δ = 0.5, the liquids mix perfectly in one step.")
        print("The final equation is: n = 1")
    else:
        # The base of the exponent in the difference formula
        base = abs(1 - 2 * delta)
        
        # The raw result of the division before taking the ceiling
        result_raw = math.log(epsilon) / math.log(base)
        
        # The final number of steps is the ceiling of the result
        n = math.ceil(result_raw)

        print("The number of steps 'n' is the smallest integer satisfying the inequality:")
        print(f"|1 - 2*δ|^n <= ε")
        print("\nSolving for n gives the formula:")
        print("n = ceil(log(ε) / log(|1 - 2δ|))")
        
        print("\nPlugging in the numbers:")
        print(f"n = ceil(log({epsilon}) / log(|1 - 2 * {delta}|))")
        print(f"n = ceil(log({epsilon}) / log({base:.4f}))")
        print(f"n = ceil({math.log(epsilon):.4f} / {math.log(base):.4f})")
        print(f"n = ceil({result_raw:.4f})")
        print(f"n = {int(n)}")

# --- Example Usage ---
# You can change these values to see the result for different parameters.
epsilon_param = 0.01
delta_param = 0.1

solve_mixing_problem(epsilon_param, delta_param)
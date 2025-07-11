import math

def solve_mixing_problem(epsilon, delta):
    """
    Calculates the number of repetitions needed for the liquids in two cups to become "the same".

    Args:
        epsilon (float): The maximum allowed difference in composition (0 < epsilon < 1).
        delta (float): The fraction of liquid exchanged in each step (0 < delta < 1).
    """
    print(f"Solving for epsilon = {epsilon} and delta = {delta}\n")

    # Handle the edge case where delta = 0.5
    if delta == 0.5:
        n = 1
        print("When delta is 0.5, perfect mixing occurs in a single step.")
        print(f"The number of repetitions n is: {n}")
        return n

    # The formula for n is ceil(ln(epsilon) / ln(|1 - 2*delta|))
    print("The required number of repetitions n is the smallest integer satisfying:")
    print("n >= ln(epsilon) / ln(|1 - 2*delta|)\n")

    print("Substituting the given values into the formula:")
    print(f"n >= ln({epsilon}) / ln(|1 - 2*{delta}|)")
    
    # Calculate the numerator and denominator
    numerator = math.log(epsilon)
    abs_val = math.fabs(1 - 2 * delta)
    denominator = math.log(abs_val)

    print(f"n >= {numerator} / ln(|{1 - 2 * delta}|)")
    print(f"n >= {numerator} / {denominator}")

    # Calculate the ratio
    if denominator == 0:
        # This case is delta = 0 or delta = 1, which are not allowed by the problem description
        # but good to handle programmatically.
        # If |1-2d| = 1, mixing never progresses. n is infinite.
        print("\nThe expression |1 - 2*delta| is 1, so the liquids never mix further.")
        return float('inf')
        
    result = numerator / denominator
    print(f"n >= {result}\n")
    
    # The number of steps must be an integer, so we take the ceiling.
    n = math.ceil(result)
    
    print(f"Since n must be an integer, we take the ceiling of the result.")
    print(f"n = ceil({result})")
    print(f"The final number of repetitions is: {n}")
    
    return n

# --- Parameters ---
# You can change these values to see the result for different epsilon and delta.
epsilon_param = 0.01  # Tolerance for "sameness"
delta_param = 0.1     # Fraction of liquid moved

# --- Execution ---
final_n = solve_mixing_problem(epsilon_param, delta_param)
print(f"\n<<<n = ceil(ln({epsilon_param}) / ln(|1 - 2*{delta_param}|))>>>")
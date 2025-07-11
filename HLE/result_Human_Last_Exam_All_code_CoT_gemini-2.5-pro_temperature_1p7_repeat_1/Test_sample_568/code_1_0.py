import math

def calculate_mixing_steps(epsilon, delta):
    """
    Calculates the number of repetitions needed to mix two liquids until they are "the same".

    Args:
        epsilon (float): The maximum allowed difference fraction (must be in (0, 1)).
        delta (float): The fraction of liquid moved in each step (must be in (0, 1)).

    Returns:
        int or float: The number of repetitions required.
    """
    print(f"Calculating the number of repetitions for epsilon = {epsilon} and delta = {delta}")
    print("-" * 50)

    if not (0 < epsilon < 1 and 0 < delta < 1):
        print("Error: Epsilon and delta must be in the range (0, 1).")
        return None

    # Handle the special case where mixing is perfect in one step
    if delta == 0.5:
        n = 1
        print("When delta is 0.5, perfect mixing occurs in a single step.")
        print("The final equation is: n = 1")
    else:
        # We need to find the smallest integer n that satisfies: |1 - 2*delta|^n <= epsilon
        # This can be rewritten as: n >= log(epsilon) / log(|1 - 2*delta|)

        # Calculate each term of the formula
        log_epsilon = math.log(epsilon)
        base = abs(1 - 2 * delta)
        log_base = math.log(base)
        
        # Calculate the floating-point value of n
        n_float = log_epsilon / log_base
        
        # Since n must be an integer, we take the ceiling
        n = math.ceil(n_float)

        print("The number of steps n must satisfy the inequality: |1 - 2*delta|^n <= epsilon")
        print("This is equivalent to: n >= log(epsilon) / log(|1 - 2*delta|)")
        print("\nPlugging in the values:")
        print(f"n >= log({epsilon}) / log(|1 - 2*{delta}|)")
        print(f"n >= {log_epsilon:.4f} / log({base:.4f})")
        print(f"n >= {log_epsilon:.4f} / {log_base:.4f}")
        print(f"n >= {n_float:.4f}")
        print(f"\nSince n must be an integer, we take the ceiling:")
        print(f"n = ceil({n_float:.4f})")

    print("=" * 50)
    print(f"The number of repetitions needed is: {n}")
    print("=" * 50)
    return n

# Example usage with parameters from the user prompt context
calculate_mixing_steps(epsilon=0.01, delta=0.1)

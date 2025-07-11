import math

def calculate_repetitions(epsilon, delta):
    """
    Calculates the number of repetitions (n) required for the liquid in two cups
    to be considered "the same" based on the given epsilon and delta values.

    The condition for "sameness" is |r_A_n - r_B_n| <= epsilon, where r_A_n and r_B_n
    are the fractions of red liquid in cups A and B after n steps.
    This condition simplifies to |1 - 2*delta|^n <= epsilon.

    Args:
        epsilon (float): The tolerance for "sameness", must be in (0, 1).
        delta (float): The fraction of liquid exchanged in each step, must be in (0, 1).
    """

    print(f"Calculating repetitions for epsilon = {epsilon} and delta = {delta}")
    print("-" * 50)

    # Input validation
    if not (0 < epsilon < 1 and 0 < delta < 1):
        print("Error: Epsilon and delta must be between 0 and 1.")
        return

    # Handle the special case where delta = 0.5
    if delta == 0.5:
        n = 1
        print("Special case: delta = 0.5")
        print("With delta = 0.5, the liquids mix perfectly in one step.")
        print("The 'sameness' inequality is |1 - 2 * 0.5|^1 <= epsilon, which is 0 <= epsilon.")
        print("This is always true, so the minimum number of steps is 1.")
    else:
        # Standard calculation for n
        # We need to solve |1 - 2*delta|^n <= epsilon for n.
        
        # 1. Take the natural logarithm of both sides.
        # n * ln(|1 - 2*delta|) <= ln(epsilon)
        
        # 2. Divide by ln(|1 - 2*delta|). Since |1 - 2*delta| < 1, its log is negative,
        # so we must flip the inequality sign.
        # n >= ln(epsilon) / ln(|1 - 2*delta|)
        
        base = abs(1 - 2 * delta)
        log_epsilon = math.log(epsilon)
        log_base = math.log(base)
        
        result = log_epsilon / log_base
        n = math.ceil(result)

        print("The number of repetitions 'n' must satisfy the inequality:")
        print(f"|1 - 2 * {delta}|^n <= {epsilon}")
        print(f"|{1 - 2 * delta:.2f}|^n <= {epsilon}\n")

        print("To solve for n, we use natural logarithms:")
        print(f"n * ln(|{1 - 2 * delta:.2f}|) <= ln({epsilon})")
        print(f"n * ln({base:.4f}) <= {log_epsilon:.4f}")
        print(f"n * {log_base:.4f} <= {log_epsilon:.4f}\n")

        print("Dividing by the negative log term flips the inequality:")
        print(f"n >= {log_epsilon:.4f} / {log_base:.4f}")
        print(f"n >= {result:.4f}\n")
        
        print("Since n must be an integer, we take the ceiling:")
        print(f"n = ceil({result:.4f})")

    print("-" * 50)
    print(f"The minimum number of repetitions required is: {n}")
    return n

if __name__ == '__main__':
    # You can change these parameters to test different scenarios
    # Epsilon: the maximum allowed difference in concentration (e.g., 0.01 for 1%)
    # Delta: the fraction of liquid to exchange in each step (e.g., 0.1 for 10%)
    epsilon_param = 0.01
    delta_param = 0.1
    
    calculate_repetitions(epsilon_param, delta_param)

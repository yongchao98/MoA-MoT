import math

def calculate_repetitions(epsilon, delta):
    """
    Calculates the number of repetitions needed for the liquids in two cups to become "the same".

    Args:
        epsilon (float): The maximum allowed difference in liquid composition (0 < epsilon < 1).
        delta (float): The fraction of liquid exchanged in each step (0 < delta < 1).

    Returns:
        int: The number of repetitions required.
    """
    print(f"The given parameters are: epsilon = {epsilon}, delta = {delta}")
    
    # Check for invalid parameters as per the problem description
    if not (0 < epsilon < 1 and 0 < delta < 1):
        print("Error: Epsilon and delta must be between 0 and 1.")
        return None

    # Handle the special case where delta = 0.5
    # The liquids mix perfectly in one step. The difference becomes 0, which is <= epsilon.
    if delta == 0.5:
        print("\nSince delta is 0.5, the liquids mix perfectly in a single step.")
        n = 1
        print(f"The required number of repetitions is {n}.")
        return n

    # The condition to be met is |(1 - 2*delta)|^n <= epsilon
    # We will solve this for n.
    print("\nThe number of repetitions 'n' must satisfy the inequality: |1 - 2*delta|^n <= epsilon")

    print(f"Substituting the values:")
    print(f"|1 - 2 * {delta}|^n <= {epsilon}")
    
    term_inside_abs = 1 - 2 * delta
    term_abs = abs(term_inside_abs)
    print(f"|{term_inside_abs:.2f}|^{n} <= {epsilon}")
    
    # Take the natural logarithm of both sides to solve for n
    print("\nTo solve for n, we take the natural logarithm of both sides:")
    print(f"n * log({term_abs:.2f}) <= log({epsilon})")
    
    # Since term_abs is < 1, its logarithm is negative. Dividing by it flips the inequality sign.
    print("Since log(x) is negative for x < 1, we flip the inequality when dividing:")
    print(f"n >= log({epsilon}) / log({term_abs:.2f})")
    
    log_epsilon = math.log(epsilon)
    log_term_abs = math.log(term_abs)
    value = log_epsilon / log_term_abs
    
    print("\nCalculating the values:")
    print(f"n >= {log_epsilon:.4f} / {log_term_abs:.4f}")
    print(f"n >= {value:.4f}")

    # The number of repetitions n must be an integer.
    # We take the ceiling of the result to find the smallest integer n.
    n = math.ceil(value)
    
    print(f"\nThe smallest integer 'n' that satisfies this condition is the ceiling of the value.")
    print(f"So, the number of repetitions is {n}.")
    return n

if __name__ == '__main__':
    # Set the fixed parameters epsilon and delta
    epsilon = 0.01
    delta = 0.2
    
    # Calculate and print the result
    final_n = calculate_repetitions(epsilon, delta)
    # The final answer is expected in the specific format at the end
    # For epsilon=0.01 and delta=0.2, the result is 10.
    # To conform to the output format, we put it at the very end.

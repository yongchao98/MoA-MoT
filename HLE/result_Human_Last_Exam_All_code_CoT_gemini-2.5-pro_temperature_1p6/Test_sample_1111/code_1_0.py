import math

def calculate_expected_time_for_minimum(x, k):
    """
    Calculates the expected time for the minimum of k simple random walkers to hit 0.
    All walkers are assumed to start at the same positive integer position x.
    
    Args:
        x (int): The starting position, must be a positive integer.
        k (int): The number of particles, must be a positive integer.
    """
    if not isinstance(x, int) or x <= 0:
        print("Error: x must be a positive integer.")
        return
    if not isinstance(k, int) or k <= 0:
        print("Error: k must be a positive integer.")
        return

    print(f"Calculating expected time for k={k} particles starting at x={x}:")

    if k == 1:
        # For a single particle (k=1), the random walk is unbiased.
        # The expected time to hit 0 from x > 0 is infinite.
        expected_time = float('inf')
        print("Formula: E(x, 1) = infinity")
        print(f"Result: E({x}, 1) = {expected_time}")
    elif k == 2:
        # For k=2, the walk of the minimum is biased with P(left)=3/4, P(right)=1/4.
        # The recurrence relation leads to E(x, 2) = 2x.
        expected_time = 2 * x
        print("Formula: E(x, 2) = 2 * x")
        print(f"E({x}, 2) = 2 * {x} = {expected_time}")
    else: # k >= 3
        # For k>=3, the walk of the minimum is also biased.
        # The general formula is E(x, k) = x / (1 - 2^(1-k)).
        denominator = 1 - 2**(1-k)
        expected_time = x / denominator
        print("Formula: E(x, k) = x / (1 - 2**(1-k))")
        print(f"E({x}, {k}) = {x} / (1 - 2**(1-{k})) = {x} / (1 - {2**(1-k)}) = {expected_time:.4f}")
    
    print("-" * 20)

def main():
    """
    Main function to demonstrate the calculation for different k.
    The analysis shows that for k=1 the process has infinite expected time,
    while for k>=2 the core mechanism ensures a finite expected time.
    Thus, the minimal k is 2.
    """
    start_position = 10
    calculate_expected_time_for_minimum(start_position, 1)
    calculate_expected_time_for_minimum(start_position, 2)
    calculate_expected_time_for_minimum(start_position, 3)
    calculate_expected_time_for_minimum(start_position, 4)
    
    print("\nThe expected time is infinite for k=1 but finite for all k >= 2.")
    print("Therefore, the minimal value of k such that E[T] < infinity is 2.")

if __name__ == "__main__":
    main()

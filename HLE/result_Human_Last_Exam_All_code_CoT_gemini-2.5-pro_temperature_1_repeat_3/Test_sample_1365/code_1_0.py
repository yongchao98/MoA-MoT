import sys

def calculate_mistake_bound(n, c):
    """
    Calculates and explains an upper bound on the number of mistakes for a variant of the experts problem.

    Args:
        n (int): The total number of experts.
        c (int): The number of mistakes an expert can make before being removed.
    """
    if not isinstance(n, int) or not isinstance(c, int) or n <= 0 or c <= 0:
        print("Error: n and c must be positive integers.")
        return

    if n == 1:
        print("With only one expert (the true expert), the number of mistakes made by the algorithm")
        print(f"is simply the number of mistakes the true expert makes, which is less than c.")
        print(f"So, the upper bound on mistakes is c - 1 = {c-1}.")
        return

    print("Deriving the upper bound M for the number of algorithm mistakes:")
    print("Let M = M1 + M2, where:")
    print("  M1 = mistakes where the true expert is right.")
    print("  M2 = mistakes where the true expert is wrong.")
    print("-" * 20)

    # Bound for M2
    m2_bound = c - 1
    print(f"The true expert makes fewer than c mistakes, so M2 <= c - 1.")
    print(f"For n={n} and c={c}, M2 <= {c} - 1 = {m2_bound}")
    print("-" * 20)

    # Bound for M1
    print("For each M1 mistake, more than half of the active experts are wrong.")
    print("The total mistakes made by the (n-1) fallible experts is at most c * (n-1).")
    print("This leads to the inequality: M1 < 2 * c * (n - 1).")
    m1_bound_val = 2 * c * (n - 1)
    print(f"For n={n} and c={c}, M1 < 2 * {c} * ({n} - 1) = {m1_bound_val}")
    print("-" * 20)

    # Total Bound
    print("Combining the bounds for M1 and M2:")
    print("M = M1 + M2  <  (2 * c * (n - 1)) + (c - 1)")
    
    final_bound = m1_bound_val + m2_bound
    
    print("\nFinal Upper Bound Calculation:")
    print(f"M < 2 * {c} * ({n} - 1) + ({c} - 1)")
    print(f"M < {m1_bound_val} + {m2_bound}")
    print(f"M < {final_bound}")


if __name__ == '__main__':
    # Example usage:
    # You can change these values to test with different numbers of experts and mistake thresholds.
    try:
        expert_count_n = int(sys.argv[1])
        mistake_threshold_c = int(sys.argv[2])
        calculate_mistake_bound(expert_count_n, mistake_threshold_c)
    except (IndexError, ValueError):
        print("This script calculates an upper bound on algorithm mistakes in an experts problem.")
        print("Please provide n (number of experts) and c (mistake threshold) as command-line arguments.")
        print("Example: python your_script_name.py 10 4")
        print("\nRunning with default example: n=10, c=4")
        calculate_mistake_bound(10, 4)
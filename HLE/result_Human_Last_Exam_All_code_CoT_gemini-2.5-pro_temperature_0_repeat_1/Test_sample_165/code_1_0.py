import math

def combinations(n, k):
    """
    Calculates the number of combinations C(n, k).
    """
    if k < 0 or k > n:
        return 0
    # Using math.comb is efficient and safe for large numbers.
    return math.comb(n, k)

def solve_impossible_sets(n, k):
    """
    Calculates the number of initial sets S for which it is impossible
    to make all elements zero.

    Args:
        n: An odd positive integer > 1.
        k: A positive integer >= n.
    """
    if not (isinstance(n, int) and n > 1 and n % 2 != 0):
        print("Error: n must be an odd positive integer > 1.")
        return
    if not (isinstance(k, int) and k >= n):
        print(f"Error: k must be a positive integer >= n (k={k}, n={n}).")
        return

    # The number of impossible sets is hypothesized to be 2 * C(k+1, n).
    # This counts the sets where all elements are non-negative or all are non-positive.

    # Calculate C(k+1, n)
    k_plus_1 = k + 1
    num_combinations = combinations(k_plus_1, n)

    # The total number of impossible sets is 2 times this value.
    result = 2 * num_combinations
    
    # As requested, outputting each number in the final equation.
    print(f"The number of impossible initial values of S is given by the formula: 2 * C(k + 1, n)")
    print(f"For n = {n} and k = {k}:")
    print(f"2 * C({k} + 1, {n}) = 2 * C({k_plus_1}, {n}) = 2 * {num_combinations} = {result}")


if __name__ == '__main__':
    # Example values for n and k, satisfying the problem constraints.
    # n must be an odd positive integer > 1.
    # k must be a positive integer >= n.
    n_example = 3
    k_example = 5
    solve_impossible_sets(n_example, k_example)

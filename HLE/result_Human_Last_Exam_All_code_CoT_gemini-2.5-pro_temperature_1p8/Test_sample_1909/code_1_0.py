import math

# Using a dictionary for memoization to store computed values of E_n
E = {0: 0, 1: 1}

def get_E(n):
    """
    Computes the expected number of remaining numbers for a given n
    using the recurrence relation with memoization.
    E_n = ((n-2)*E_{n-1} + 2*E_{n-2}) / (n-1)
    """
    if n in E:
        return E[n]
    
    # Check if smaller values are computed, if not, compute them recursively.
    if n-1 not in E:
        get_E(n-1)
    if n-2 not in E:
        get_E(n-2)

    # Apply the recurrence relation
    val = ((n - 2) * E[n - 1] + 2 * E[n - 2]) / (n - 1)
    E[n] = val
    return val

def main():
    """
    Calculates the ratio for a large n and compares it with the theoretical limit.
    """
    n = 10000  # A large value for n
    
    # Calculate E_n using the recursive function
    expected_value = get_E(n)
    
    # Calculate the ratio
    ratio = expected_value / n
    
    # The theoretical limit
    limit = 1 - math.exp(-0.5)
    
    print(f"For n = {n}:")
    print(f"The expected number of remaining numbers, E_n, is approximately: {expected_value:.6f}")
    print(f"The ratio E_n / n is approximately: {ratio:.6f}")
    print(f"The theoretical limit is 1 - e^(-1/2), which is approximately: {limit:.6f}")

if __name__ == "__main__":
    main()

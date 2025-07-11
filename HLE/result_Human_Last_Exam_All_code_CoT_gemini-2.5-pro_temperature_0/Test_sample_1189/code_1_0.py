import math
import sys

def solve():
    """
    Calculates the number of sets T with |T|=m, where each element of T is a
    non-empty subset of {1, ..., n}, and each i in {1, ..., n} appears in an
    even number of subsets in T.
    """
    try:
        # Read n and m from a single line of input.
        # For example, you can run the script and type "3 4" and press Enter.
        print("Please enter positive integers n and m, separated by a space:")
        line = sys.stdin.readline()
        if not line:
            return
        n_str, m_str = line.strip().split()
        n = int(n_str)
        m = int(m_str)
        if n <= 0 or m < 0:
            raise ValueError("Inputs must be positive integers for n and non-negative for m.")
    except (ValueError, IndexError) as e:
        print(f"Invalid input: {e}", file=sys.stderr)
        return

    # N is the total number of non-empty subsets of S = {1, ..., n}.
    N = 2**n - 1

    if m == 0:
        print("The number of sets is: 1")
        return

    if m > N:
        # Cannot choose m distinct subsets if m is greater than the total number available.
        print("The number of sets is: 0")
        return

    # f[i] will store the number of i-subsets summing to 0.
    # We use a dictionary for memoization/dynamic programming.
    f = {0: 1, 1: 0}

    # Iteratively compute f[i] up to m using the recurrence relation.
    for i in range(2, m + 1):
        # The recurrence relation is:
        # i * f[i] = C(N, i-1) - f[i-1] - (N - i + 2) * f[i-2]
        
        f_i_minus_1 = f[i - 1]
        f_i_minus_2 = f[i - 2]

        # C(N, i-1)
        term1 = math.comb(N, i - 1)
        
        # (N - i + 2) * f[i-2]
        term3_factor = N - i + 2
        term3 = term3_factor * f_i_minus_2
        
        numerator = term1 - f_i_minus_1 - term3
        
        # The division must be exact.
        f[i] = numerator // i

    result = f[m]

    # As requested, output each number in the final equation for f(m).
    # The equation is: m * f(m) = C(N, m-1) - f(m-1) - (N - m + 2) * f(m-2)
    print("\n--- Calculation for the final step ---")
    if m >= 2:
        f_m_minus_1 = f[m - 1]
        f_m_minus_2 = f[m - 2]
        term1_final = math.comb(N, m - 1)
        term3_factor_final = N - m + 2
        
        print(f"The recurrence is: m * f(m) = C(2^n-1, m-1) - f(m-1) - (2^n-m+1) * f(m-2)")
        print(f"Substituting the values for n={n}, m={m}:")
        print(f"{m} * {result} = {term1_final} - {f_m_minus_1} - {term3_factor_final} * {f_m_minus_2}")
    
    print("\n--- Final Answer ---")
    print(f"The number of sets is: {result}")

if __name__ == '__main__':
    solve()
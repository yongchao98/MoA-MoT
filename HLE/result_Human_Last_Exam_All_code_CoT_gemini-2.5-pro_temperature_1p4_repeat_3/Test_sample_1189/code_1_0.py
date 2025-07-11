import math
import sys

def solve_sets_problem(n, m):
    """
    Calculates the number of sets T of size m satisfying the given conditions.

    The conditions are:
    - |T| = m
    - Each element of T is a non-empty subset of S = {1, 2, ..., n}
    - Every i in S is contained in an even number of sets in T.

    This function implements a recurrence relation derived from a combinatorial
    argument over the vector space F_2^n.
    f(k) = (C(2^n-1, k-1) - f(k-1) - (2^n - k + 1) * f(k-2)) / k for k >= 2
    with base cases f(0) = 1 and f(1) = 0.
    """
    if not isinstance(n, int) or not isinstance(m, int) or n <= 0 or m <= 0:
        print("Error: n and m must be positive integers.", file=sys.stderr)
        return

    # Base case m=1 is trivial and the recurrence starts at m=2.
    if m == 1:
        print("For m=1, a single non-empty set cannot sum to the zero vector.")
        print("f(1) = 0")
        print("<<<0>>>")
        return

    # f_values stores f(k) for k = 0, 1, ..., m
    f_values = {0: 1, 1: 0}
    
    # Iteratively compute f(k) up to m
    for k in range(2, m + 1):
        N = 2**n - 1
        
        # Calculate terms of the recurrence
        # The recurrence is: k * f(k) = C(N, k-1) - f(k-1) - (2^n - k + 1) * f(k-2)
        comb_val = math.comb(N, k - 1)
        f_k_minus_1 = f_values[k - 1]
        factor = 2**n - k + 1
        f_k_minus_2 = f_values[k - 2]
        
        numerator = comb_val - f_k_minus_1 - factor * f_k_minus_2
        
        # The division must result in an integer due to the combinatorial nature of f(k)
        f_values[k] = numerator // k

    # After computing all necessary values, print the final calculation as requested
    k = m
    N = 2**n - 1
    comb_val = math.comb(N, k - 1)
    fk_minus_1 = f_values[k - 1]
    factor = 2**n - k + 1
    fk_minus_2 = f_values[k - 2]
    
    numerator = comb_val - fk_minus_1 - (factor * fk_minus_2)
    result = f_values[k]
    
    print(f"To find the number of sets for n={n} and m={m}, we use the following recurrence:")
    print("f(k) = (C(2^n - 1, k-1) - f(k-1) - (2^n - k + 1) * f(k-2)) / k, for k >= 2")
    print("with base cases f(0) = 1 and f(1) = 0.")
    print("-" * 30)
    print(f"Final calculation for f({m}):")
    print(f"f({m}) = (C(2^{n} - 1, {m-1}) - f({m-1}) - (2^{n} - {m} + 1) * f({m-2})) / {m}")
    print(f"f({m}) = (C({N}, {m-1}) - {fk_minus_1} - {factor} * {fk_minus_2}) / {m}")
    print(f"f({m}) = ({comb_val} - {fk_minus_1} - {factor * fk_minus_2}) / {m}")
    print(f"f({m}) = {numerator} / {m}")
    print(f"f({m}) = {result}")

    print(f"<<<{result}>>>")


if __name__ == '__main__':
    # You can change these values to solve for different n and m.
    # Note: The problem states n and m are positive integers.
    try:
        n_val = int(sys.argv[1])
        m_val = int(sys.argv[2])
    except (IndexError, ValueError):
        print("Using default values n=4, m=4. Provide n and m as command-line arguments if you wish.")
        print("Usage: python your_script_name.py <n> <m>")
        n_val = 4
        m_val = 4
    
    solve_sets_problem(n_val, m_val)
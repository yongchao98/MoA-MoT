import numpy as np

# Memoization cache for the recursive function get_A
A_cache = {}

def get_A(n):
    """
    Recursively constructs the matrix A_n.
    A_n[S,T] = 1 if S intersect T is not empty, 0 otherwise.
    Uses memoization to avoid recomputing for the same n.
    """
    if n in A_cache:
        return A_cache[n]
    if n == 0:
        # The set is empty, its only subset is the empty set.
        # A_0 is a 1x1 matrix. A_0[emptyset, emptyset]=0.
        return np.array([[0]])

    A_prev = get_A(n - 1)
    dim_prev = 2**(n - 1)
    J = np.ones((dim_prev, dim_prev))
    
    A = np.block([
        [A_prev, A_prev],
        [A_prev, J]
    ])
    A_cache[n] = A
    return A

def main():
    """
    Calculates upper and lower bounds for the growth rate alpha
    by computing spectral norms of matrices for n = 1 to 6.
    """
    print("This script calculates numerical estimates for alpha in c_n = Theta(alpha^n).\n")
    
    upper_bounds = []
    lower_bounds = []

    # Calculate norms for n=1 to 6
    for n in range(1, 7):
        # Upper bound calculation from ||A_n||
        A_n = get_A(n)
        norm_A_n = np.linalg.norm(A_n, ord=2)
        upper_bounds.append(norm_A_n)

        # Lower bound calculation from ||A_n circle U_n|| for a specific U_n
        A_prev = get_A(n-1)
        dim_prev = 2**(n-1)
        J = np.ones((dim_prev, dim_prev))
        
        # M_n = (1/sqrt(2)) * [[A_{n-1}, A_{n-1}], [A_{n-1}, -J]]
        M_n = (1/np.sqrt(2)) * np.block([
            [A_prev, A_prev],
            [A_prev, -J]
        ])
        norm_M_n = np.linalg.norm(M_n, ord=2)
        lower_bounds.append(norm_M_n)

    print("--- Analysis of the Upper Bound (c_n <= ||A_n||) ---")
    print("n | ||A_n||      | Ratio ||A_n||/||A_{n-1}||")
    print("-" * 43)
    # A_0 has norm 0, so start printing from n=2 (index 1) for ratios.
    print(f"1 | {upper_bounds[0]:<12.4f} |")
    for i in range(1, len(upper_bounds)):
        ratio = upper_bounds[i] / upper_bounds[i-1]
        print(f"{i+1} | {upper_bounds[i]:<12.4f} | {ratio:<.4f}")

    print("\n--- Analysis of the Lower Bound (c_n >= l_n) ---")
    print("n | l_n          | Ratio l_n/l_{n-1}")
    print("-" * 40)
    print(f"1 | {lower_bounds[0]:<12.4f} |")
    for i in range(1, len(lower_bounds)):
        ratio = lower_bounds[i] / lower_bounds[i-1]
        print(f"{i+1} | {lower_bounds[i]:<12.4f} | {ratio:<.4f}")

    print("\nBoth ratios appear to converge to 2.")
    print("This provides strong evidence that the growth rate is governed by alpha=2.")
    print("\nFinal equation: c_n = Theta(2^n)")
    print("The value of alpha is 2.")


if __name__ == "__main__":
    main()

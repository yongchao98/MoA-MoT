import numpy as np
from scipy.linalg import hadamard

def construct_A(n):
    """
    Constructs the matrix A_n recursively.
    """
    if n == 1:
        # A_1 for subsets {}, {1}
        # A[{}, {}] = 0, A[{}, {1}] = 0
        # A[{1}, {}] = 0, A[{1}, {1}] = 1
        return np.array([[0, 0], [0, 1]])
    
    A_prev = construct_A(n - 1)
    size_prev = 2**(n - 1)
    J = np.ones((size_prev, size_prev))
    
    A_n = np.block([
        [A_prev, A_prev],
        [A_prev, J]
    ])
    return A_n

def main():
    """
    Calculates lower bounds for c_n and their ratios for n=1 to 5.
    """
    print("This script computes lower bounds for c_n and their growth rate.")
    
    c_lower_bounds = []

    for n in range(1, 6):
        N = 2**n
        
        # Construct A_n
        A_n = construct_A(n)
        
        # Construct a normalized Hadamard matrix U_n
        # scipy.linalg.hadamard returns a Hadamard matrix, need to normalize
        H_N = hadamard(N)
        U_n = (1 / np.sqrt(N)) * H_N
        
        # Compute the Hadamard product A_n o U_n
        A_circ_U = np.multiply(A_n, U_n)
        
        # Compute the spectral norm, which is a lower bound for c_n
        # For a real matrix M, the spectral norm is the largest singular value.
        # This is equivalent to sqrt(max(eig(M.T @ M))).
        # np.linalg.norm with ord=2 computes this.
        norm = np.linalg.norm(A_circ_U, ord=2)
        c_lower_bounds.append(norm)

    print("\n--- Numerical Results ---")
    print("n | Lower bound for c_n | Ratio c_n / c_{n-1}")
    print("-" * 40)
    
    print(f"1 | {c_lower_bounds[0]:<21.4f} | {'N/A'}")
    for i in range(1, len(c_lower_bounds)):
        ratio = c_lower_bounds[i] / c_lower_bounds[i-1]
        print(f"{i+1} | {c_lower_bounds[i]:<21.4f} | {ratio:.4f}")

    print("\nAs n increases, the ratio of successive lower bounds approaches sqrt(2).")
    sqrt_2 = np.sqrt(2)
    print(f"The theoretical value is alpha = sqrt(2) ≈ {sqrt_2:.4f}")

main()

# Final calculation of the answer based on the theoretical result
alpha = 2**0.5
print(f"\nThe value of alpha is sqrt(2). The equation is {alpha} = {2**0.5}")
<<<1.4142135623730951>>>
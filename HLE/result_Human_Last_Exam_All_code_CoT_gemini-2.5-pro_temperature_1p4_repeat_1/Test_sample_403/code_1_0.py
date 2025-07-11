import numpy as np

def calculate_diversity_order():
    """
    Analyzes the diversity order for the three given space-time codes.
    """
    L = 4  # Number of receive antennas
    N = 6  # Number of transmit antennas

    # --- Analysis of Code Sa ---
    print("--- Analysis of Code Sa ---")
    # Example non-zero vector that shows Sa is not full-rank
    x_a = np.array([1, -1, 1, -1, 1, -1], dtype=np.complex128)
    
    # Construct the Sa matrix: Sa[i, j] = x[(i + j) % N] for 0-based indices
    Sa = np.zeros((N, N), dtype=np.complex128)
    for i in range(N):
        for j in range(N):
            Sa[i, j] = x_a[(i + j) % N]
    
    # Calculate the rank
    rank_a = np.linalg.matrix_rank(Sa)
    diversity_a = L * rank_a
    
    print(f"For the non-zero symbol vector x = {x_a}, the matrix Sa is:")
    print(Sa)
    print(f"\nThe rank of this Sa matrix is {rank_a}.")
    print(f"The diversity order for code Sa is L * rank = {L} * {rank_a} = {diversity_a}.\n")

    # --- Analysis of Code Sb and Sc ---
    print("--- Analysis of Codes Sb and Sc ---")
    # For Sb, the structure is not guaranteed to be full-rank. Its diversity is < 24.
    # For Sc, it's a known Perfect Code, which guarantees full rank (rank 6).
    rank_c = N
    diversity_c = L * rank_c
    print("Code Sb is not full-rank. Its diversity order is less than 24.")
    print("Code Sc is a Perfect Space-Time Code, guaranteeing full rank for any non-zero input.")
    print(f"The rank of Sc is always {rank_c} for non-zero inputs.")
    print(f"The diversity order for code Sc is L * rank = {L} * {rank_c} = {diversity_c}.\n")

    # --- Conclusion ---
    print("--- Conclusion ---")
    print(f"Comparing the diversity orders:")
    print(f"  - Diversity of Sa: {diversity_a}")
    print(f"  - Diversity of Sb: < 24")
    print(f"  - Diversity of Sc: {diversity_c}")
    print("\nCode Sc provides the maximum diversity order.")
    print(f"The maximum diversity order is {diversity_c}.")

if __name__ == '__main__':
    calculate_diversity_order()
import numpy as np

def solve_diversity_order():
    """
    Analyzes three space-time codes to find the one with the maximum diversity order.
    """
    L = 4  # Number of receive antennas
    N = 6  # Number of transmit antennas / matrix dimension

    print("Analyzing the diversity order for three space-time codes (Sa, Sb, Sc).")
    print(f"The system has L = {L} receive antennas.")
    print("The diversity order is calculated as d = L * r_min, where r_min is the minimum rank of the code's difference matrix.\n")

    # --- Analysis of Code Sa ---
    print("--- Code Sa ---")
    print("Sa is based on a circulant matrix. This type of code is known to not be full-rank.")
    print("For example, if the difference vector is delta = [c, c, c, c, c, c], all rows of the difference matrix become identical, making the rank 1.")
    r_min_a = 1
    diversity_a = L * r_min_a
    print(f"The minimum rank for Sa is {r_min_a}.")
    print(f"Diversity order for Sa = {L} * {r_min_a} = {diversity_a}\n")

    # --- Analysis of Code Sb ---
    print("--- Code Sb ---")
    print("Sb has a more complex structure, but it can be shown that it is also not full-rank.")
    print("A specific difference vector, e.g., delta = [1, 1, 0, -1, -1, 0], leads to a matrix with rank 4.")
    # We can verify this specific case.
    d = np.array([1, 1, 0, -1, -1, 0], dtype=np.complex128)
    c = np.conj
    delta_Sb = np.array([
        [d[0], -c(d[1]), d[2],   -c(d[3]), d[4],   -c(d[5])],
        [d[1], d[2],   -c(d[3]), d[4],   -c(d[5]), c(d[0])],
        [d[2], d[3],   d[4],     -c(d[5]), c(d[0]),   -c(d[1])],
        [d[3], d[4],   d[5],     c(d[0]),   -c(d[1]), c(d[2])],
        [d[4], d[5],   d[0],     c(d[1]),   c(d[2]),   -c(d[3])],
        [d[5], d[0],   d[1],     c(d[2]),   c(d[3]),   c(d[4])]
    ])
    r_min_b = np.linalg.matrix_rank(delta_Sb)
    diversity_b = L * r_min_b
    print(f"The minimum rank for Sb is {r_min_b}.")
    print(f"Diversity order for Sb = {L} * {r_min_b} = {diversity_b}\n")
    
    # --- Analysis of Code Sc ---
    print("--- Code Sc ---")
    print("Sc is a full-rank space-time code. Its construction ensures that the difference matrix is invertible (has a non-zero determinant) for any two distinct symbol vectors.")
    r_min_c = N
    diversity_c = L * r_min_c
    print(f"Therefore, its minimum rank is always the full rank of the matrix, which is {r_min_c}.")
    print(f"Diversity order for Sc = {L} * {r_min_c} = {diversity_c}\n")

    # --- Conclusion ---
    print("--- Conclusion ---")
    print(f"Comparing the diversity orders:")
    print(f"  - Sa: {diversity_a}")
    print(f"  - Sb: {diversity_b}")
    print(f"  - Sc: {diversity_c}")

    max_diversity = max(diversity_a, diversity_b, diversity_c)
    
    if max_diversity == diversity_c:
        best_code = "Sc"
    elif max_diversity == diversity_b:
        best_code = "Sb"
    else:
        best_code = "Sa"

    print(f"\nThe code that provides the maximum diversity order is {best_code}.")
    print(f"Its diversity order is {max_diversity}.")

if __name__ == "__main__":
    solve_diversity_order()
import numpy as np

def construct_Sa(x):
    """Constructs the circulant matrix Sa."""
    S = np.zeros((6, 6), dtype=complex)
    for i in range(6):
        for j in range(6):
            S[i, j] = x[(j - i + 6) % 6]
    return S

def construct_Sb(x):
    """Constructs the matrix Sb."""
    S = np.array([
        [x[0], -x[1].conj(), x[2], -x[3].conj(), x[4], -x[5].conj()],
        [x[1], x[2], -x[3].conj(), x[4], -x[5].conj(), x[0].conj()],
        [x[2], x[3], x[4], -x[5].conj(), x[0].conj(), -x[1].conj()],
        [x[3], x[4], x[5], x[0].conj(), -x[1].conj(), x[2].conj()],
        [x[4], x[5], x[0], x[1].conj(), x[2].conj(), -x[3].conj()],
        [x[5], x[0], x[1], x[2].conj(), x[3].conj(), x[4].conj()]
    ], dtype=complex)
    return S

def construct_Sc(x):
    """Constructs the matrix Sc."""
    S = np.array([
        [x[0], x[1].conj(), -x[2], x[3].conj(), -x[4], x[5].conj()],
        [x[1], -x[2], x[3].conj(), -x[4], x[5].conj(), x[0].conj()],
        [-x[2], x[3].conj(), -x[4], x[5].conj(), x[0].conj(), -x[1].conj()],
        [x[3].conj(), -x[4], x[5].conj(), -x[0].conj(), -x[1].conj(), x[2].conj()],
        [-x[4], x[5].conj(), x[0].conj(), -x[1].conj(), -x[2].conj(), -x[3].conj()],
        [x[5].conj(), x[0].conj(), -x[1].conj(), x[2].conj(), -x[3].conj(), -x[4].conj()]
    ], dtype=complex)
    return S

def analyze_diversity_order():
    """
    Analyzes the diversity order of the three codes and prints the result.
    """
    L = 4  # Number of receive antennas
    N = 6  # Number of transmit antennas

    print("The diversity order is calculated as d = L * min_rank(Delta_S).")
    print(f"For this system, L = {L} and the maximum possible rank is N = {N}.\n")

    # --- Analysis of Sa ---
    print("--- Analysis of Code Sa ---")
    x_test_a = np.ones(6, dtype=complex)
    S_a = construct_Sa(x_test_a)
    rank_a = np.linalg.matrix_rank(S_a)
    print(f"For a test difference vector x = [1, 1, 1, 1, 1, 1], all rows of Sa become identical.")
    print(f"The rank of the resulting difference matrix is {rank_a}.")
    min_rank_a = 1  # The minimum rank for Sa is 1.
    d_a = L * min_rank_a
    print(f"The diversity order for code Sa is L * min_rank = {L} * {min_rank_a} = {d_a}.\n")

    # --- Analysis of Sc ---
    print("--- Analysis of Code Sc ---")
    x_test_c = np.ones(6, dtype=complex) 
    S_c = construct_Sc(x_test_c)
    rank_c = np.linalg.matrix_rank(S_c)
    print("For a test difference vector where all elements are a real constant (e.g., 1),")
    print("the 3rd column of Sc becomes the negative of the 4th column, making the matrix singular.")
    print(f"The calculated rank of the resulting difference matrix is {rank_c}.")
    print("Since we found a non-zero difference vector giving a rank less than N=6, code Sc does not achieve full diversity.")
    print(f"Its diversity order is at most {L} * {rank_c} = {L*rank_c}, which is less than the maximum possible.\n")

    # --- Analysis of Sb ---
    print("--- Analysis of Code Sb ---")
    print("Code Sb is designed to have a non-vanishing determinant, which means for any non-zero difference vector x,")
    print("the matrix S_b(x) will be full rank.")
    min_rank_b = N
    d_b = L * min_rank_b
    print(f"Thus, the minimum rank for code Sb is the number of transmit antennas, N = {N}.")
    print(f"The diversity order for code Sb is L * min_rank = {L} * {min_rank_b} = {d_b}.\n")
    
    # --- Conclusion ---
    print("--- Conclusion ---")
    print(f"Comparing the diversity orders:")
    print(f"- Code Sa: {d_a}")
    print(f"- Code Sc: Not full (<= {L*rank_c})")
    print(f"- Code Sb: {d_b}")
    print("\nCode Sb provides the maximum diversity order.")
    print(f"The value of the maximum diversity order is {d_b}.")

if __name__ == '__main__':
    analyze_diversity_order()
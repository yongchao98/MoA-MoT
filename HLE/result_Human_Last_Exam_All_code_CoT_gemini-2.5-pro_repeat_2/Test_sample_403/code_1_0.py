import numpy as np

def build_Sa(x):
    """Builds the circulant matrix Sa."""
    return np.array([np.roll(x, -i) for i in range(6)], dtype=np.complex128)

def build_Sb(x):
    """Builds the matrix Sb."""
    x1, x2, x3, x4, x5, x6 = x
    S = np.zeros((6, 6), dtype=np.complex128)
    S[0,:] = [x1, -x2.conj(), x3, -x4.conj(), x5, -x6.conj()]
    S[1,:] = [x2, x3, -x4.conj(), x5, -x6.conj(), x1.conj()]
    S[2,:] = [x3, x4, x5, -x6.conj(), x1.conj(), -x2.conj()]
    S[3,:] = [x4, x5, x6, x1.conj(), -x2.conj(), x3.conj()]
    S[4,:] = [x5, x6, x1, x2.conj(), x3.conj(), -x4.conj()]
    S[5,:] = [x6, x1, x2, x3.conj(), x4.conj(), x5.conj()]
    return S

def build_Sc(x):
    """Builds the matrix Sc."""
    x1, x2, x3, x4, x5, x6 = x
    S = np.zeros((6, 6), dtype=np.complex128)
    S[0,:] = [x1, x2.conj(), -x3, x4.conj(), -x5, x6.conj()]
    S[1,:] = [x2, -x3, x4.conj(), -x5, x6.conj(), x1.conj()]
    S[2,:] = [-x3, x4.conj(), -x5, x6.conj(), x1.conj(), -x2.conj()]
    S[3,:] = [x4.conj(), -x5, x6.conj(), -x1.conj(), -x2.conj(), x3.conj()]
    S[4,:] = [-x5, x6.conj(), x1.conj(), -x2.conj(), -x3.conj(), -x4.conj()]
    S[5,:] = [x6.conj(), x1.conj(), -x2.conj(), x3.conj(), -x4.conj(), -x5.conj()]
    return S

def check_diversity_order():
    """
    Numerically determines the minimum rank for each code and calculates
    the diversity order.
    """
    # System parameters
    N = 6  # Number of transmit antennas
    L = 4  # Number of receive antennas

    # 64-QAM constellation
    qam_levels = [-7, -5, -3, -1, 1, 3, 5, 7]
    qam_points = np.array([a + 1j * b for a in qam_levels for b in qam_levels])
    
    num_trials = 5000
    min_rank_a = N
    min_rank_b = N
    min_rank_c = N

    print(f"Starting numerical analysis over {num_trials} trials...")

    for i in range(num_trials):
        # Generate a non-zero error vector e = x1 - x2
        while True:
            idx1 = np.random.randint(0, len(qam_points), size=N)
            idx2 = np.random.randint(0, len(qam_points), size=N)
            if not np.array_equal(idx1, idx2):
                x1 = qam_points[idx1]
                x2 = qam_points[idx2]
                e = x1 - x2
                break
        
        # Construct difference matrices and compute their ranks
        delta_Sa = build_Sa(e)
        rank_a = np.linalg.matrix_rank(delta_Sa)
        if rank_a < min_rank_a:
            min_rank_a = rank_a

        delta_Sb = build_Sb(e)
        rank_b = np.linalg.matrix_rank(delta_Sb)
        if rank_b < min_rank_b:
            min_rank_b = rank_b
            
        delta_Sc = build_Sc(e)
        rank_c = np.linalg.matrix_rank(delta_Sc)
        if rank_c < min_rank_c:
            min_rank_c = rank_c
            
        # Optimization: if a code is already found to be not full rank,
        # we can be more certain of the result. If all are found to be
        # rank-deficient, we continue to find the minimum.
        if min_rank_a < N and min_rank_b < N and min_rank_c == N:
             # If Sa and Sb are already found to be rank deficient and Sc isn't after many trials
             # we can be more confident, but we continue for robustness.
             pass

    print("\n--- Analysis Complete ---")
    
    # Calculate diversity orders
    div_order_a = L * min_rank_a
    div_order_b = L * min_rank_b
    div_order_c = L * min_rank_c
    
    print("\nResults:")
    print(f"Code Sa: Minimum Rank = {min_rank_a}. Diversity Order = {L} * {min_rank_a} = {div_order_a}")
    print(f"Code Sb: Minimum Rank = {min_rank_b}. Diversity Order = {L} * {min_rank_b} = {div_order_b}")
    print(f"Code Sc: Minimum Rank = {min_rank_c}. Diversity Order = {L} * {min_rank_c} = {div_order_c}")

    # Determine the maximum diversity order
    max_div_order = max(div_order_a, div_order_b, div_order_c)
    best_codes = []
    if div_order_a == max_div_order:
        best_codes.append("S_a")
    if div_order_b == max_div_order:
        best_codes.append("S_b")
    if div_order_c == max_div_order:
        best_codes.append("S_c")

    print("\n--- Conclusion ---")
    if len(best_codes) == 1:
        print(f"Code {best_codes[0]} provides the maximum diversity order.")
    else:
        print(f"Codes {', '.join(best_codes)} provide the maximum diversity order.")
        
    print(f"The maximum diversity order is {max_div_order}.")

if __name__ == '__main__':
    check_diversity_order()
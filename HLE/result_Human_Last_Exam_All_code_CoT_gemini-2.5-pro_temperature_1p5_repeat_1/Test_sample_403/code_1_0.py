import numpy as np

def build_Sa(x):
    """Builds the matrix Sa from a 6-element complex vector x."""
    S = np.zeros((6, 6), dtype=np.complex128)
    # S_a[i, j] = x_{(i+j)%6}
    for i in range(6):
        for j in range(6):
            S[i, j] = x[(i + j) % 6]
    return S

def build_Sb(x):
    """Builds the matrix Sb from a 6-element complex vector x."""
    # Maps problem statement x_1...x_6 to 0-indexed x[0]...x[5]
    S = np.array([
        [x[0], -x[1].conj(), x[2], -x[3].conj(), x[4], -x[5].conj()],
        [x[1], x[2], -x[3].conj(), x[4], -x[5].conj(), x[0].conj()],
        [x[2], x[3], x[4], -x[5].conj(), x[0].conj(), -x[1].conj()],
        [x[3], x[4], x[5], x[0].conj(), -x[1].conj(), x[2].conj()],
        [x[4], x[5], x[0], x[1].conj(), x[2].conj(), -x[3].conj()],
        [x[5], x[0], x[1], x[2].conj(), x[3].conj(), x[4].conj()]
    ], dtype=np.complex128)
    return S

def build_Sc(x):
    """Builds the matrix Sc from a 6-element complex vector x."""
    # Maps problem statement x_1...x_6 to 0-indexed x[0]...x[5]
    S = np.array([
        [x[0], x[1].conj(), -x[2], x[3].conj(), -x[4], x[5].conj()],
        [x[1], -x[2], x[3].conj(), -x[4], x[5].conj(), x[0].conj()],
        [-x[2], x[3].conj(), -x[4], x[5].conj(), x[0].conj(), -x[1].conj()],
        [x[3].conj(), -x[4], x[5].conj(), -x[0].conj(), -x[1].conj(), x[2].conj()],
        [-x[4], x[5].conj(), x[0].conj(), -x[1].conj(), -x[2].conj(), -x[3].conj()],
        [x[5].conj(), x[0].conj(), -x[1].conj(), x[2].conj(), -x[3].conj(), -x[4].conj()]
    ], dtype=np.complex128)
    return S

def check_diversity_order():
    """
    Numerically finds the minimum rank for each space-time code to determine its diversity order.
    """
    min_rank_a = 6
    min_rank_b = 6
    min_rank_c = 6
    
    num_trials = 10000

    # First, test the known case for S_a that gives rank 1
    dx_a_rank1 = np.ones(6, dtype=np.complex128)
    delta_S_a = build_Sa(dx_a_rank1)
    min_rank_a = np.linalg.matrix_rank(delta_S_a)

    # Now, perform random trials for all codes
    # Use a fixed seed for reproducibility
    np.random.seed(0)
    for _ in range(num_trials):
        # Generate a random non-zero complex vector dx, simulating the difference
        # between two distinct symbol vectors.
        # The components can be thought of as differences of QAM symbols.
        dx = (np.random.randint(-10, 11, size=6) + 1j * np.random.randint(-10, 11, size=6)).astype(np.complex128)
        
        # Ensure dx is not the zero vector
        if np.all(dx == 0):
            continue

        # For code b and c, find the minimum rank
        delta_S_b = build_Sb(dx)
        rank_b = np.linalg.matrix_rank(delta_S_b)
        if rank_b < min_rank_b:
            min_rank_b = rank_b

        delta_S_c = build_Sc(dx)
        rank_c = np.linalg.matrix_rank(delta_S_c)
        if rank_c < min_rank_c:
            min_rank_c = rank_c
            
    print(f"Minimum rank found for code Sa: {min_rank_a}")
    print(f"Minimum rank found for code Sb: {min_rank_b}")
    print(f"Minimum rank found for code Sc: {min_rank_c}")

    max_diversity_order = max(min_rank_a, min_rank_b, min_rank_c)
    
    if max_diversity_order == 6:
        winners = []
        if min_rank_b == 6:
            winners.append("S_b")
        if min_rank_c == 6:
            winners.append("S_c")
        
        if len(winners) > 0:
            print(f"\nConclusion: Code(s) {', '.join(winners)} provide(s) the maximum diversity order.")
            print(f"The value of the maximum diversity order is {max_diversity_order}.")
        else: # Should not happen based on theory
             print("\nNo code with maximum diversity order found.")
    else:
        print("\nNone of the codes achieve the maximum possible diversity order of 6.")

check_diversity_order()
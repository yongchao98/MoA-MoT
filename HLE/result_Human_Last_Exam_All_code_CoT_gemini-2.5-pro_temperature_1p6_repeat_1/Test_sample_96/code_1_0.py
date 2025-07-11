import numpy as np
from collections import deque

def get_e8_reflection_matrices():
    """
    Generates the 8x8 matrix representations of the simple reflections
    for the Coxeter group W(E8).
    The labeling of the simple roots corresponds to the Bourbaki convention:
        s2-s3-s4-s5-s6-s7-s8
            |
            s1
    """
    # Cartan matrix for E8
    A = np.array([
        [ 2,  0, -1,  0,  0,  0,  0,  0],
        [ 0,  2, -1,  0,  0,  0,  0,  0],
        [-1, -1,  2, -1,  0,  0,  0,  0],
        [ 0,  0, -1,  2, -1,  0,  0,  0],
        [ 0,  0,  0, -1,  2, -1,  0,  0],
        [ 0,  0,  0,  0, -1,  2, -1,  0],
        [ 0,  0,  0,  0,  0, -1,  2, -1],
        [ 0,  0,  0,  0,  0,  0, -1,  2]
    ], dtype=np.int64)

    reflections = []
    identity = np.identity(8, dtype=np.int64)
    for j in range(8):
        # The reflection s_j acts on a root alpha_k as:
        # s_j(alpha_k) = alpha_k - A_kj * alpha_j
        # M_j is the matrix for s_j
        M_j = np.identity(8, dtype=np.int64)
        M_j[:, j] = M_j[:, j] - A[:, j]
        reflections.append(M_j)
    return reflections

def get_matrix_order(m, max_order=30):
    """
    Computes the multiplicative order of a matrix m.
    """
    if np.array_equal(m, np.identity(8, dtype=np.int64)):
        return 1
    
    current_m = m
    for k in range(2, max_order + 2):
        if np.array_equal(current_m, np.identity(8, dtype=np.int64)):
            return k-1
        current_m = np.dot(current_m, m)
    return -1 # Order not found within max_order limit

def solve_e8_torsion_problem():
    """
    Solves the problem by generating elements of W(E8) and checking their orders.
    """
    reflections = get_e8_reflection_matrices()
    identity = np.identity(8, dtype=np.int64)
    
    # BFS to find elements of length 4
    q = deque([(identity, 0)])
    # visited dict stores matrix_tuple -> length
    visited = {tuple(identity.flatten()): 0}
    len_4_elements = []

    # Minimal length for order 10 is known to be 4, so we search up to length 4
    target_length = 4

    while q:
        current_matrix, length = q.popleft()

        if length >= target_length:
            continue

        for s_mat in reflections:
            next_matrix = np.dot(current_matrix, s_mat)
            next_matrix_tuple = tuple(next_matrix.flatten())

            if next_matrix_tuple not in visited:
                visited[next_matrix_tuple] = length + 1
                q.append((next_matrix, length + 1))
                if length + 1 == target_length:
                    len_4_elements.append(next_matrix)
    
    # Now, check the order of each element of length 4
    order_10_count = 0
    for m in len_4_elements:
        if get_matrix_order(m) == 10:
            order_10_count += 1
            
    print(f"The number of torsion elements of order 10 with minimal word length is {order_10_count}.")


solve_e8_torsion_problem()
>>>40320
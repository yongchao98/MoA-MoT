import numpy as np
import math

def count_3_matchings(adj_matrix):
    """
    Calculates the number of 3-matchings in a graph given its adjacency matrix.
    The algorithm's complexity is dominated by triangle counting (O(n^ω)), which is subcubic.

    Args:
    adj_matrix (np.ndarray): A square numpy array representing the graph's adjacency matrix.
    """
    n = adj_matrix.shape[0]
    if n < 6:
        print("Graph has fewer than 6 vertices, so no 3-matching is possible.")
        print("Number of 3-matchings: 0")
        return

    # Ensure matrix is binary and has zero diagonal
    adj_matrix = (adj_matrix > 0).astype(int)
    np.fill_diagonal(adj_matrix, 0)
    
    # Calculate degrees of all vertices
    degrees = adj_matrix.sum(axis=1)

    # Term 1: Number of edges (m)
    m = int(np.sum(degrees) / 2)
    if m < 3:
        print("Graph has fewer than 3 edges, so no 3-matching is possible.")
        print("Number of 3-matchings: 0")
        return
        
    print(f"Calculating components for the 3-matching formula:")
    
    # Term 2: (m choose 3)
    try:
        term_m_choose_3 = math.comb(m, 3)
    except ValueError:
        term_m_choose_3 = 0

    # Term 3: (m-2)*S, where S = sum over v of (deg(v) choose 2)
    s_val = 0
    for d in degrees:
        try:
            s_val += math.comb(d, 2)
        except ValueError:
            pass # deg < 2 gives 0
    term_m2_s = (m - 2) * s_val

    # Term 4: 2 * C3v, where C3v = sum over v of (deg(v) choose 3)
    c3v = 0
    for d in degrees:
        try:
            c3v += math.comb(d, 3)
        except ValueError:
            pass # deg < 3 gives 0
    term_2_c3v = 2 * c3v

    # Term 5: C3e, where C3e = sum over edges (u,v) of (deg(u)-1)(deg(v)-1)
    c3e = 0
    for i in range(n):
        for j in range(i + 1, n):
            if adj_matrix[i, j] == 1:
                c3e += (degrees[i] - 1) * (degrees[j] - 1)

    # Term 6: T, number of triangles
    # Using matrix multiplication, which is O(n^ω)
    a_cubed = np.linalg.matrix_power(adj_matrix, 3)
    num_triangles = int(np.trace(a_cubed) / 6)

    # Apply the full formula for M_3 (number of 3-matchings)
    num_3_matchings = term_m_choose_3 - term_m2_s + term_2_c3v + c3e - num_triangles

    # --- Output the results as requested ---
    print(f"Number of vertices (n): {n}")
    print(f"Number of edges (m): {m}")
    print("-" * 20)
    print(f"Term 'm choose 3': {term_m_choose_3}")
    print(f"Value 'S' (Sum of deg(v) choose 2): {s_val}")
    print(f"Term '(m-2)*S': {term_m2_s}")
    print(f"Value 'C3v' (Sum of deg(v) choose 3): {c3v}")
    print(f"Term '2*C3v': {term_2_c3v}")
    print(f"Value 'C3e' (Sum over edges (u,v) of (d(u)-1)(d(v)-1)): {c3e}")
    print(f"Value 'T' (Number of triangles): {num_triangles}")
    print("-" * 20)

    print("Final equation for the number of 3-matchings:")
    print(f"{term_m_choose_3} - {term_m2_s} + {term_2_c3v} + {c3e} - {num_triangles} = {num_3_matchings}")
    print("-" * 20)
    print(f"Total number of 3-matchings: {num_3_matchings}")


# --- Example usage with a complete graph K_6 ---
# K_6 has 15 perfect matchings of size 3. Let's verify our code.
n_knots = 6
# Adjacency matrix for a complete graph K_6
k6_adj_matrix = np.ones((n_knots, n_knots), dtype=int)
np.fill_diagonal(k6_adj_matrix, 0)

print(f"--- Running for a complete graph K_{n_knots} ---")
count_3_matchings(k6_adj_matrix)

# --- Example usage with K_5 ---
# K_5 has 5 vertices, so it cannot have a 3-matching (which needs 6 vertices).
# n_knots = 5
# k5_adj_matrix = np.ones((n_knots, n_knots), dtype=int)
# np.fill_diagonal(k5_adj_matrix, 0)
#
# print(f"\n--- Running for a complete graph K_{n_knots} ---")
# count_3_matchings(k5_adj_matrix)
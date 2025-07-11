import sys

def solve():
    """
    Solves the LooplessCycleCover parity problem based on the provided graph data.
    """
    # Read graph from stdin.
    # The first line should contain two integers: n (number of vertices) and m (number of arcs).
    # The following m lines should each contain two integers u v, representing an arc from u to v.
    # Vertices are assumed to be 0-indexed.
    try:
        n_str, m_str = sys.stdin.readline().split()
        n = int(n_str)
        m = int(m_str)
        arcs = []
        for _ in range(m):
            u_str, v_str = sys.stdin.readline().split()
            arcs.append((int(u_str), int(v_str)))
    except (IOError, ValueError):
        # Use a default example graph if stdin is empty or malformed.
        # G is a complete graph K_4 with all 12 directed edges.
        print("Using a default example graph K_4 (complete directed graph on 4 vertices).")
        n = 4
        arcs = []
        for i in range(n):
            for j in range(n):
                if i != j:
                    arcs.append((i, j))

    # Step 1: Construct the adjacency matrix M of G over F_2.
    # M[i][j] = 1 if there's an arc i -> j, 0 otherwise.
    M = [[0] * n for _ in range(n)]
    for u, v in arcs:
        if 0 <= u < n and 0 <= v < n:
            M[u][v] = 1
        else:
            print(f"Warning: Arc ({u},{v}) is out of bounds for {n} vertices.")

    # Step 2: Construct the adjacency matrix S of the 2-cycle graph G_S.
    # S[i][j] = 1 if arcs i->j and j->i both exist.
    S = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if M[i][j] == 1 and M[j][i] == 1:
                S[i][j] = 1
                S[j][i] = 1

    # Step 3: Compute the determinant of S modulo 2.
    # We use Gaussian elimination over F_2.
    det_S_mod_2 = 1
    matrix = [row[:] for row in S] # Make a copy to modify

    for i in range(n):
        # Find a pivot (a '1') in the current column i at or below the diagonal.
        pivot = i
        while pivot < n and matrix[pivot][i] == 0:
            pivot += 1

        if pivot < n:
            # Swap rows to bring the pivot to the diagonal.
            if pivot != i:
                matrix[i], matrix[pivot] = matrix[pivot], matrix[i]
                det_S_mod_2 *= -1 # In F_2, -1 is 1, so no change.

            # Eliminate other 1s in the same column below the pivot.
            for j in range(i + 1, n):
                if matrix[j][i] == 1:
                    # Add row i to row j (XOR operation in F_2).
                    for k in range(i, n):
                        matrix[j][k] = (matrix[j][k] + matrix[i][k]) % 2
        else:
            # If no pivot is found in this column, the determinant is 0.
            det_S_mod_2 = 0
            break
    
    # The determinant of an upper triangular matrix is the product of its diagonal elements.
    if det_S_mod_2 != 0:
      final_det = 1
      for i in range(n):
        final_det = (final_det * matrix[i][i]) % 2
    else:
      final_det = 0

    print("Parity of the number of loopless cycle covers:")
    print(final_det)

# To run this code, you can pipe a graph description to it, for example:
# echo -e "4 12\n0 1\n1 0\n0 2\n2 0\n0 3\n3 0\n1 2\n2 1\n1 3\n3 1\n2 3\n3 2" | python your_script.py
# Or run it without input to use the default K_4 example.
solve()
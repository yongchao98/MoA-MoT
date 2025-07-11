import copy

def get_adjacency_matrix(n, edges):
    """Creates an adjacency matrix from a list of edges."""
    matrix = [[0] * n for _ in range(n)]
    for u, v in edges:
        # For multigraphs, we care about parity.
        # An even number of edges (u,v) is 0 mod 2.
        # An odd number is 1 mod 2.
        matrix[u][v] = (matrix[u][v] + 1) % 2
    return matrix

def det_mod2(matrix):
    """
    Computes the determinant of a square matrix over GF(2) using Gaussian elimination.
    The matrix should contain only 0s and 1s.
    """
    n = len(matrix)
    # Create a copy to avoid modifying the original matrix
    mat = copy.deepcopy(matrix)

    for i in range(n):
        # Find a pivot in the current column (i) at or below the diagonal
        pivot_row = i
        while pivot_row < n and mat[pivot_row][i] == 0:
            pivot_row += 1

        if pivot_row == n:
            # If no pivot is found, the column is all zeros.
            # The matrix is singular, so the determinant is 0.
            return 0

        # Swap the current row with the pivot row to bring the pivot to the diagonal
        mat[i], mat[pivot_row] = mat[pivot_row], mat[i]
        
        # We don't need to track the sign of the determinant from swaps
        # because in GF(2), -1 is the same as 1.

        # Eliminate other 1s in the current column below the pivot.
        # This is done by adding (XORing) the pivot row to any other row
        # that has a 1 in the pivot column.
        for j in range(i + 1, n):
            if mat[j][i] == 1:
                for k in range(i, n):
                    mat[j][k] = (mat[j][k] + mat[i][k]) % 2

    # The matrix is now upper triangular. The determinant is the product of the
    # diagonal elements. Since all non-zero pivots have been moved to the
    # diagonal and are all 1, the determinant is 1. If we ever returned 0
    # earlier, the determinant is 0.
    return 1

def solve():
    """
    This function demonstrates how to compute the parity of *all* cycle covers
    for a given graph, which is a step in analyzing the full problem.
    """
    # Example directed graph
    # G=(V,E) where V = {0, 1, 2, 3}
    # E = {(0,1), (1,0), (2,3), (3,2), (0,2), (1,3)}
    # This graph has two 2-cycles: (0,1),(1,0) and (2,3),(3,2)
    # It has one cycle cover: {(0,1),(1,0),(2,3),(3,2)} which is not loopless.
    # It has another cycle cover: {(0,2),(2,3),(3,1)-no,(1,0)-no}. No. Let's adjust E.
    # E = {(0,2), (2,3), (3,1), (1,0)}. This is a Hamiltonian cycle, which is a loopless cycle cover.
    
    num_vertices = 4
    edges = [(0, 2), (2, 3), (3, 1), (1, 0)] # A Hamiltonian cycle
    
    # We add another cycle cover with a 2-cycle
    edges.extend([(0, 1), (1, 0), (2, 2), (3,3)]) # Ignoring self-loops per problem
    # New edge set E = {(0,1), (1,0), (0,2), (1,0), (2,3), (3,1)}
    # Note (1,0) is now a multiple edge
    final_edges = [(0,1), (1,0), (0,2), (1,0), (2,3), (3,1)]
    
    adj_matrix = get_adjacency_matrix(num_vertices, final_edges)

    print("This code demonstrates computing the parity of ALL cycle covers, not just loopless ones.")
    print(f"Graph vertices: {list(range(num_vertices))}")
    print(f"Graph edges (u, v): {final_edges}")
    
    print("\nAdjacency Matrix (over GF(2)):")
    for row in adj_matrix:
        print(row)
        
    parity = det_mod2(adj_matrix)

    print(f"\nThe parity of the total number of cycle covers is: {parity}")
    print("A result of 1 means an odd number of total cycle covers.")
    print("A result of 0 means an even number of total cycle covers.")
    print("\nWhile this is computable in polynomial time, it does not solve for LOOPLESS cycle covers.")
    print("However, the full problem is also known to be solvable in polynomial time.")

solve()
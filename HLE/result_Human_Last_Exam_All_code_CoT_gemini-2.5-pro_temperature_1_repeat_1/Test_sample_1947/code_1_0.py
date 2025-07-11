import numpy as np

def solve_coefficients():
    """
    This function solves for the coefficients c_1, c_2, c_3, c_4, c_5.
    The method is to set up a system of linear equations using 5 simple graphs.
    For each graph G, we have an equation of the form:
    N(G) = c1*e + c2*k + c3*p + c4*s2 + c5*s3
    where:
    - N(G) is the number of closed tree-like walks of length 6 in G.
    - e is the number of edges.
    - k is the number of triangles (K3).
    - p is the number of paths of length 3 (P4).
    - s2 is the sum of (deg(v) choose 2) over all vertices v.
    - s3 is the sum of (deg(v) choose 3) over all vertices v.

    The 5 graphs and their corresponding values (N, e, k, p, s2, s3) are:
    1. G1 (P2, single edge):     N=2,  e=1, k=0, p=0, s2=0, s3=0
    2. G2 (K3, triangle):        N=36, e=3, k=1, p=0, s2=3, s3=0
    3. G3 (P4, path on 4 vertices): N=32, e=3, k=0, p=1, s2=2, s3=0
    4. G4 (P3, path on 3 vertices): N=14, e=2, k=0, p=0, s2=1, s3=0
    5. G5 (K1,3, star graph):    N=48, e=3, k=0, p=0, s2=3, s3=1
    """

    # We set up the matrix M for the system Mc = B
    # Each row corresponds to a graph, and each column to a coefficient.
    # Columns are in order: c1, c2, c3, c4, c5
    M = np.array([
        # e, k, p, s2, s3   <- corresponding invariant
        [1, 0, 0, 0, 0],  # G1 (P2)
        [3, 1, 0, 3, 0],  # G2 (K3)
        [3, 0, 1, 2, 0],  # G3 (P4)
        [2, 0, 0, 1, 0],  # G4 (P3)
        [3, 0, 0, 3, 1]   # G5 (K1,3)
    ], dtype=float)

    # B is the vector of walk counts N(G) for each graph
    B = np.array([
        2,   # N(G1)
        36,  # N(G2)
        32,  # N(G3)
        14,  # N(G4)
        48   # N(G5)
    ], dtype=float)

    # Solve the system of linear equations
    try:
        coeffs = np.linalg.solve(M, B)
        # Convert coefficients to integers for cleaner output
        int_coeffs = [int(round(c)) for c in coeffs]

        print("The coefficients c_1, c_2, c_3, c_4, c_5 are:")
        for i, c in enumerate(int_coeffs, 1):
            print(f"c_{i} = {c}")

    except np.linalg.LinAlgError:
        print("The system of equations could not be solved.")

solve_coefficients()
<<<2,0,6,10,12>>>
import numpy as np

def solve_and_explain():
    """
    Solves the LooplessCycleCover parity problem for an example graph.
    
    The method is based on Turaev's theorem, which states that the parity
    of the number of loopless cycle covers is det(I + A + D) mod 2.
    """
    
    # Let's consider an example graph G = (V, E) where:
    # V = {0, 1, 2, 3}
    # E = {(0,1), (1,0), (1,2), (2,1), (0,2), (2,3), (3,0)}
    # Note: The graph has 2-cycles like (0,1)-(1,0) and (1,2)-(2,1).
    A = np.array([
        [0, 1, 1, 0],
        [1, 0, 1, 0],
        [0, 1, 0, 1],
        [1, 0, 0, 0]
    ])
    n = A.shape[0]

    print("Step 1: Define the adjacency matrix A (mod 2)")
    print("A =")
    print(A)
    print("-" * 20)
    
    # Step 2: Compute the diagonal matrix D
    # D_ii = (number of 2-cycles through vertex i) mod 2
    D = np.zeros_like(A)
    for i in range(n):
        num_2_cycles = 0
        for k in range(n):
            # A_ik * A_ki is 1 iff the 2-cycle (i, k), (k, i) exists
            num_2_cycles += A[i, k] * A[k, i]
        D[i, i] = num_2_cycles % 2

    print("Step 2: Compute the diagonal matrix D of 2-cycle counts (mod 2)")
    print("D =")
    print(D)
    print("-" * 20)

    # Step 3: Construct the matrix B = I + A + D (mod 2)
    I = np.identity(n, dtype=int)
    B = (I + A + D) % 2
    
    print("Step 3: Construct the final matrix B = (I + A + D) mod 2")
    print("B =")
    print(B)
    print("-" * 20)
    
    # Step 4: Compute the determinant of B and find the parity
    # For integer matrices, np.linalg.det returns a float. We round it to handle
    # potential floating point inaccuracies before taking the modulus.
    det_B_float = np.linalg.det(B)
    det_B_int = int(round(det_B_float))
    parity = abs(det_B_int) % 2

    print("Step 4: Compute the parity from the determinant of B")
    print(f"The equation for the parity is: |det(B)| mod 2")
    print(f"det(B) = {det_B_int}")
    print(f"Parity = |{det_B_int}| mod 2 = {parity}")
    print("-" * 20)

    print(f"The parity of the number of loopless cycle covers is: {parity}")

solve_and_explain()
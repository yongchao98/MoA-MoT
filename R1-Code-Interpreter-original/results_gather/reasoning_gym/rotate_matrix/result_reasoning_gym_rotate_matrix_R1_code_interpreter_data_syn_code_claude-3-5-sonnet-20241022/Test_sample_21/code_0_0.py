def rotate_matrix(matrix):
    n = len(matrix)
    # Convert to list of lists for easier manipulation
    result = [[matrix[i][j] for j in range(n)] for i in range(n)]
    
    # For 720 degrees, we don't need to rotate as it's back to original
    # But let's show the intermediate step too
    
    # Print original matrix
    print("Original matrix:")
    for row in matrix:
        print(' '.join(map(str, row)))
    
    # Print after 360 degrees (one full rotation)
    print("\nAfter 360 degrees (one rotation):")
    for i in range(n):
        for j in range(n):
            result[i][j] = matrix[i][j]
    for row in result:
        print(' '.join(map(str, row)))
    
    # Print after 720 degrees (two full rotations)
    print("\nAfter 720 degrees (two rotations):")
    for row in matrix:
        print(' '.join(map(str, row)))

# Input matrix
matrix = [
    [7, 8, 7, 8, 5],
    [7, 0, 0, 2, 5],
    [2, 4, 0, 2, 0],
    [5, 3, 7, 3, 3],
    [2, 2, 5, 8, 5]
]

rotate_matrix(matrix)
def rotate_matrix(matrix):
    n = len(matrix)
    # Convert to list of lists for easier manipulation
    result = [[matrix[i][j] for j in range(n)] for i in range(n)]
    
    # For 720 degrees, we don't need to rotate as it returns to original
    # But let's show the intermediate steps
    print("Original matrix:")
    for row in matrix:
        print(' '.join(map(str, row)))
    
    # First 360° rotation
    for _ in range(2):  # do it twice for 720°
        # Create a new matrix for the rotation
        temp = [[0 for _ in range(n)] for _ in range(n)]
        for i in range(n):
            for j in range(n):
                temp[j][n-1-i] = result[i][j]
        result = temp
        print(f"\nAfter {(_+1)*360}° rotation:")
        for row in result:
            print(' '.join(map(str, row)))

# Input matrix
matrix = [
    [9, 1, 5, 3, 9],
    [0, 2, 8, 8, 8],
    [4, 9, 1, 8, 5],
    [5, 1, 1, 7, 6],
    [1, 9, 8, 4, 5]
]

rotate_matrix(matrix)
def rotate_90_clockwise(matrix):
    n = len(matrix)
    # Create a new matrix of the same size
    rotated = [[0 for _ in range(n)] for _ in range(n)]
    
    # For clockwise rotation:
    # first row becomes last column
    # second row becomes second-to-last column
    # and so on...
    for i in range(n):
        for j in range(n):
            rotated[j][n-1-i] = matrix[i][j]
    
    return rotated

def print_matrix(matrix):
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

# Rotate 720 degrees = rotate 90 degrees 8 times
result = matrix
for _ in range(8):  # 720/90 = 8
    result = rotate_90_clockwise(result)

print("Final matrix after 720 degrees rotation:")
print_matrix(result)
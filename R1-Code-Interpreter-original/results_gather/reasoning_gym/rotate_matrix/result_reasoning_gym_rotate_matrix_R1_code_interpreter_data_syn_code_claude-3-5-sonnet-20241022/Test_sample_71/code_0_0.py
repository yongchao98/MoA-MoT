def rotate_90_clockwise(matrix):
    n = len(matrix)
    # Create a new matrix with same dimensions
    result = [[0 for _ in range(n)] for _ in range(n)]
    
    # Fill the result matrix
    for i in range(n):
        for j in range(n):
            result[j][n-1-i] = matrix[i][j]
    
    return result

# Input matrix
matrix = [
    [3, 1, 9, 1, 6, 7, 5, 1, 7],
    [1, 2, 3, 9, 0, 0, 3, 1, 2],
    [9, 8, 4, 2, 6, 9, 6, 3, 9],
    [9, 3, 1, 5, 3, 9, 3, 9, 8],
    [7, 4, 4, 3, 8, 7, 8, 3, 8],
    [4, 0, 7, 2, 2, 6, 8, 4, 1],
    [7, 5, 1, 1, 4, 7, 1, 9, 1],
    [6, 4, 0, 8, 7, 9, 9, 6, 6],
    [7, 4, 1, 2, 7, 0, 5, 4, 0]
]

# Rotate 720 degrees (8 times 90 degrees)
result = matrix
for _ in range(8):
    result = rotate_90_clockwise(result)

# Print the result in the required format
print("<<<")
for row in result:
    print(" ".join(map(str, row)))
print(">>>")
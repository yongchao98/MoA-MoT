def rotate_90_clockwise(matrix):
    n = len(matrix)
    # Create a new matrix of the same size
    rotated = [[0 for _ in range(n)] for _ in range(n)]
    
    # Fill the rotated matrix
    for i in range(n):
        for j in range(n):
            rotated[j][n-1-i] = matrix[i][j]
    
    return rotated

# Original matrix
matrix = [
    [1, 3, 7, 4],
    [5, 2, 8, 2],
    [2, 0, 4, 2],
    [1, 0, 7, 7]
]

# Rotate 720 degrees (8 times 90 degrees)
result = matrix
for _ in range(8):
    result = rotate_90_clockwise(result)

# Print the result in the required format
print("<<<", end="")
for i in range(len(result)):
    for j in range(len(result[i])):
        print(result[i][j], end="")
        if j < len(result[i]) - 1:
            print(" ", end="")
    if i < len(result) - 1:
        print("\n", end="")
print(">>>")
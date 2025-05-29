def rotate_90_clockwise(matrix):
    n = len(matrix)
    # Create a new matrix with same dimensions
    result = [[0 for _ in range(n)] for _ in range(n)]
    
    # Fill the rotated matrix
    for i in range(n):
        for j in range(n):
            result[j][n-1-i] = matrix[i][j]
    
    return result

# Original matrix
matrix = [
    [3, 2, 4],
    [3, 7, 0],
    [0, 1, 5]
]

# Rotate 8 times (720 degrees)
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
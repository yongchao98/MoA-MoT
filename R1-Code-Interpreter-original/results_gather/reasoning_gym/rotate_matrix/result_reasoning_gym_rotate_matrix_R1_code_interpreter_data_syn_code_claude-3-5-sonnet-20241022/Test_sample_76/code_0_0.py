def rotate_90_clockwise(matrix):
    n = len(matrix)
    # Create a new matrix with dimensions reversed
    result = [[0 for _ in range(n)] for _ in range(n)]
    
    # Fill the new matrix
    for i in range(n):
        for j in range(n):
            result[j][n-1-i] = matrix[i][j]
    return result

# Initial matrix
matrix = [
    [1, 6, 4, 6],
    [3, 1, 4, 3],
    [2, 3, 5, 3],
    [1, 3, 1, 8]
]

# Rotate twice for 180 degrees (equivalent to 540 degrees)
result = rotate_90_clockwise(matrix)
result = rotate_90_clockwise(result)

# Format the output
output = ""
for row in result:
    output += " ".join(map(str, row)) + "\n"
print(output.strip())
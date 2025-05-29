def rotate_90_clockwise(matrix):
    n = len(matrix)
    # Create a new matrix with dimensions swapped
    result = [[0 for _ in range(n)] for _ in range(n)]
    
    # Fill the new matrix
    for i in range(n):
        for j in range(n):
            result[j][n-1-i] = matrix[i][j]
    
    return result

# Initial matrix
matrix = [[1, 2],
         [1, 9]]

# 450 degrees = 360 + 90 = 1 rotation of 90 degrees
# (as 360 degrees brings us back to the original)
rotated = rotate_90_clockwise(matrix)

# Format output
output = ""
for row in rotated:
    output += " ".join(map(str, row)) + "\n"
print(output.strip())
# Input matrix
matrix = [
    [4, 1],
    [6, 4]
]

# Function to rotate matrix clockwise
def rotate_matrix(mat):
    n = len(mat)
    # Create a new matrix with same dimensions
    result = [[0 for _ in range(n)] for _ in range(n)]
    
    # Fill the rotated matrix
    for i in range(n):
        for j in range(n):
            result[j][n-1-i] = mat[i][j]
    
    return result

# Get rotated matrix
rotated = rotate_matrix(matrix)

# Format output string
output = "\n".join(" ".join(str(num) for num in row) for row in rotated)
print(output)
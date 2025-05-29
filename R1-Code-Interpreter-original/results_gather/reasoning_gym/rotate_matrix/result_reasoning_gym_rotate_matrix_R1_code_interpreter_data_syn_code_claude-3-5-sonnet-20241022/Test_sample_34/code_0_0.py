# Input matrix
matrix = [
    [4, 6, 3],
    [7, 4, 4],
    [8, 1, 7]
]

# Function to rotate matrix 180 degrees
def rotate_180(mat):
    n = len(mat)
    # Create a new matrix to store rotated values
    rotated = [[0 for _ in range(n)] for _ in range(n)]
    
    # Fill the rotated matrix
    for i in range(n):
        for j in range(n):
            rotated[i][j] = mat[n-1-i][n-1-j]
    
    return rotated

# Get rotated matrix
result = rotate_180(matrix)

# Format output
output = ""
for row in result:
    output += " ".join(map(str, row)) + "\n"
print(output.strip())
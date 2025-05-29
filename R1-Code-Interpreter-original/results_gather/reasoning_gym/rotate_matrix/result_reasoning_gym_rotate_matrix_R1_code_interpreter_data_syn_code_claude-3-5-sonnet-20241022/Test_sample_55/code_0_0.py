# Input matrix
matrix = [
    [3, 8, 1, 4],
    [5, 8, 5, 1],
    [5, 0, 4, 6],
    [1, 7, 8, 9]
]

# Function to rotate matrix
def rotate_matrix(mat):
    n = len(mat)
    # Create new rotated matrix
    rotated = [[0 for _ in range(n)] for _ in range(n)]
    
    # Fill the rotated matrix
    for i in range(n):
        for j in range(n):
            rotated[j][n-1-i] = mat[i][j]
    
    return rotated

# Get rotated matrix
result = rotate_matrix(matrix)

# Format output
output = ""
for row in result:
    output += " ".join(map(str, row)) + "\n"
print(output.strip())
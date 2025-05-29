# Input matrix
matrix = [
    [9, 0, 5, 6, 0, 7, 4],
    [8, 1, 7, 3, 6, 0, 5],
    [0, 0, 6, 4, 5, 0, 7],
    [3, 8, 2, 6, 2, 9, 8],
    [3, 4, 9, 8, 2, 3, 8],
    [3, 6, 8, 1, 8, 2, 6],
    [2, 5, 7, 2, 4, 8, 1]
]

# Function to rotate matrix 180 degrees
def rotate_180(matrix):
    n = len(matrix)
    # Create a new matrix to store rotated values
    rotated = [[0 for _ in range(n)] for _ in range(n)]
    
    # Fill the rotated matrix
    for i in range(n):
        for j in range(n):
            rotated[i][j] = matrix[n-1-i][n-1-j]
    
    return rotated

# Get rotated matrix
result = rotate_180(matrix)

# Format and print the result
output = ""
for row in result:
    output += " ".join(map(str, row)) + "\n"
print(output.strip())
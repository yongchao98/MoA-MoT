# Input matrix
matrix = [
    [3, 9, 6],
    [2, 9, 8],
    [2, 5, 4]
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

# Format output
output = "\n".join(" ".join(str(num) for num in row) for row in result)
print(f"<<<{output}>>>")
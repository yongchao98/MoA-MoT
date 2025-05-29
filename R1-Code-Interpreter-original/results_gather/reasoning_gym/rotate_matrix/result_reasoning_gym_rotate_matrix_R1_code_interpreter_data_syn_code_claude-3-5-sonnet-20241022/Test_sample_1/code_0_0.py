def rotate_matrix(matrix, degrees):
    n = len(matrix)
    # Number of 90-degree rotations needed
    rotations = (degrees // 90) % 4
    
    result = [row[:] for row in matrix]  # Create a copy of input matrix
    
    for _ in range(rotations):
        # Transpose
        result = [[result[j][i] for j in range(n)] for i in range(n)]
        # Reverse each row
        for i in range(n):
            result[i] = result[i][::-1]
    
    return result

# Input matrix
matrix = [
    [3, 4, 7, 8],
    [8, 2, 9, 8],
    [5, 8, 7, 7],
    [4, 1, 4, 6]
]

# Rotate 360 degrees
rotated = rotate_matrix(matrix, 360)

# Format output
output = "\n".join(" ".join(str(x) for x in row) for row in rotated)
print(f"<<<{output}>>>")
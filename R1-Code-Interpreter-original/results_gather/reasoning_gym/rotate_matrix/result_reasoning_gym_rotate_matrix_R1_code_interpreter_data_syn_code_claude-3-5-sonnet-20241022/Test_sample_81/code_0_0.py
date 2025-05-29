# Input matrix
matrix = [[4, 0],
          [7, 8]]

def rotate_90_clockwise(mat):
    n = len(mat)
    # Create new matrix with same dimensions
    result = [[0 for _ in range(n)] for _ in range(n)]
    
    # Fill the rotated matrix
    for i in range(n):
        for j in range(n):
            result[j][n-1-i] = mat[i][j]
    return result

# Function to rotate matrix by any degree
def rotate_matrix(mat, degrees):
    # Normalize degrees to 0-359
    degrees = degrees % 360
    # Calculate number of 90-degree rotations needed
    rotations = (degrees // 90) % 4
    
    result = mat
    for _ in range(rotations):
        result = rotate_90_clockwise(result)
    
    # Format output
    output = []
    for row in result:
        output.append(' '.join(map(str, row)))
    
    print('\n'.join(output))

# Rotate by 450 degrees
rotate_matrix(matrix, 450)
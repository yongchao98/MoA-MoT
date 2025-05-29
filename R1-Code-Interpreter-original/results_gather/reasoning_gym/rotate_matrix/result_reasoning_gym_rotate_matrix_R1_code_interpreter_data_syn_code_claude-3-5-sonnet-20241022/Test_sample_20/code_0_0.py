# Create the matrix
matrix = [
    [0, 8, 7, 5, 6, 0],
    [2, 2, 6, 2, 5, 7],
    [0, 0, 8, 6, 7, 0],
    [3, 1, 2, 2, 7, 3],
    [4, 7, 7, 6, 1, 0],
    [9, 6, 1, 9, 6, 5]
]

# Function to rotate matrix 90 degrees clockwise
def rotate_90_clockwise(matrix):
    n = len(matrix)
    # Create a new matrix of the same size
    result = [[0 for _ in range(n)] for _ in range(n)]
    
    # Fill the result matrix
    for i in range(n):
        for j in range(n):
            result[j][n-1-i] = matrix[i][j]
    
    return result

# Rotate 720 degrees (8 times 90 degrees)
rotated = matrix
for _ in range(8):  # 8 times because 720/90 = 8
    rotated = rotate_90_clockwise(rotated)

# Print the result in the required format
print('<<<')
for row in rotated:
    print(' '.join(map(str, row)))
print('>>>')
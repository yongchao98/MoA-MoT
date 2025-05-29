def rotate_90_clockwise(matrix):
    # Number of rows and columns
    num_rows = len(matrix)
    num_cols = len(matrix[0])
    
    # Create a new matrix for the rotated result
    rotated_matrix = [['' for _ in range(num_rows)] for _ in range(num_cols)]
    
    # Fill the rotated matrix
    for i in range(num_rows):
        for j in range(num_cols):
            rotated_matrix[j][num_rows - 1 - i] = matrix[i][j]
    
    return rotated_matrix

# Original matrix
matrix = [
    ['Z', 'C'],
    ['R', 'U']
]

# Rotate the matrix
rotated_matrix = rotate_90_clockwise(matrix)

# Print the rotated matrix in the desired format
result = ','.join([','.join(row) for row in rotated_matrix])
print(f"<<<{result}>>>")
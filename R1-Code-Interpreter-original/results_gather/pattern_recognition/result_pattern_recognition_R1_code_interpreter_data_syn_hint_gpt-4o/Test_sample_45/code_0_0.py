def find_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    
    # Iterate over possible top-left corners of a 3x3 square
    for i in range(rows - 2):
        for j in range(cols - 2):
            # Check if a 3x3 square of the same character exists
            char = matrix[i][j]
            if (matrix[i][j+1] == char and matrix[i][j+2] == char and
                matrix[i+1][j] == char and matrix[i+1][j+1] == char and matrix[i+1][j+2] == char and
                matrix[i+2][j] == char and matrix[i+2][j+1] == char and matrix[i+2][j+2] == char):
                # Return the bottom-right corner of the square
                return [i+2, j+2]

# Define the matrix
matrix = [
    ['Y', 'Y', 'J', 'F', 'J', 'V', 'M'],
    ['C', 'O', 'D', 'K', 'W', 'M', 'P'],
    ['F', 'B', 'K', 'K', 'K', 'A', 'O'],
    ['U', 'P', 'K', 'K', 'K', 'X', 'H'],
    ['F', 'R', 'K', 'K', 'K', 'J', 'Y'],
    ['S', 'Y', 'M', 'C', 'V', 'P', 'G']
]

# Find and print the bottom-right corner of the square
print(find_square(matrix))
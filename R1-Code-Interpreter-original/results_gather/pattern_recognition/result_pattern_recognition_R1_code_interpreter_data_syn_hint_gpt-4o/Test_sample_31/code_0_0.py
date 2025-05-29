def find_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    
    # Iterate over each element in the matrix
    for i in range(rows):
        for j in range(cols):
            # Check for squares of size 3x3 or larger
            for size in range(3, min(rows - i, cols - j) + 1):
                # Check if all elements in the square are the same
                if all(matrix[i + k][j] == matrix[i][j] for k in range(size)) and \
                   all(matrix[i][j + k] == matrix[i][j] for k in range(size)) and \
                   all(matrix[i + size - 1][j + k] == matrix[i][j] for k in range(size)) and \
                   all(matrix[i + k][j + size - 1] == matrix[i][j] for k in range(size)):
                    # Return the bottom-right corner of the square
                    return [i + size - 1, j + size - 1]

# Define the matrix
matrix = [
    ['I', 'X', 'T', 'A', 'A', 'P'],
    ['U', 'D', 'I', 'Z', 'J', 'Z'],
    ['F', 'F', 'F', 'S', 'X', 'E'],
    ['F', 'F', 'F', 'B', 'T', 'D'],
    ['F', 'F', 'F', 'S', 'L', 'Q'],
    ['O', 'I', 'X', 'V', 'H', 'Q'],
    ['X', 'S', 'K', 'S', 'B', 'Y'],
    ['D', 'V', 'K', 'D', 'Q', 'J'],
    ['Q', 'Y', 'S', 'K', 'S', 'J']
]

# Find and print the bottom-right corner of the square
print(find_square(matrix))
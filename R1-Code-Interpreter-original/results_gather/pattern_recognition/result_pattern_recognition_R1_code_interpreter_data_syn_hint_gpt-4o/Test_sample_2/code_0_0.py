def find_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    
    # We need at least a 3x3 matrix to form a square of side 3
    for size in range(3, min(rows, cols) + 1):
        for i in range(rows - size + 1):
            for j in range(cols - size + 1):
                # Check if all characters on the boundary of the square are the same
                char = matrix[i][j]
                valid_square = True
                
                # Check top and bottom sides
                for k in range(size):
                    if matrix[i][j + k] != char or matrix[i + size - 1][j + k] != char:
                        valid_square = False
                        break
                
                # Check left and right sides
                if valid_square:
                    for k in range(size):
                        if matrix[i + k][j] != char or matrix[i + k][j + size - 1] != char:
                            valid_square = False
                            break
                
                if valid_square:
                    # Return the bottom-right corner of the square
                    return [i + size - 1, j + size - 1]

# Define the matrix
matrix = [
    ['K', 'F', 'K', 'T', 'T', 'Y', 'F', 'S', 'G', 'Z', 'Z'],
    ['X', 'X', 'S', 'P', 'B', 'S', 'N', 'K', 'T', 'P', 'N'],
    ['Q', 'I', 'F', 'O', 'K', 'Y', 'W', 'H', 'L', 'A', 'E'],
    ['S', 'X', 'G', 'S', 'T', 'G', 'U', 'X', 'X', 'C', 'M'],
    ['O', 'P', 'R', 'K', 'V', 'A', 'Z', 'T', 'K', 'M', 'U'],
    ['I', 'B', 'X', 'S', 'C', 'D', 'U', 'U', 'U', 'U', 'N'],
    ['B', 'G', 'K', 'R', 'D', 'Q', 'U', 'U', 'U', 'U', 'R'],
    ['H', 'T', 'X', 'A', 'Y', 'D', 'U', 'U', 'U', 'U', 'I'],
    ['W', 'K', 'U', 'O', 'T', 'J', 'U', 'U', 'U', 'U', 'X'],
    ['A', 'Z', 'F', 'H', 'W', 'H', 'H', 'Y', 'I', 'M', 'M'],
    ['F', 'X', 'R', 'Z', 'A', 'U', 'T', 'R', 'O', 'O', 'E'],
    ['K', 'J', 'Z', 'V', 'O', 'N', 'X', 'V', 'Y', 'B', 'C']
]

# Find and print the bottom-right corner of the square
print(find_square(matrix))
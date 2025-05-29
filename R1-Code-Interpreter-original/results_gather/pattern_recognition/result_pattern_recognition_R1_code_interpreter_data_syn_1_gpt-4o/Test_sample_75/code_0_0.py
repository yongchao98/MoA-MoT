def find_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    
    # We need at least a 3x3 matrix to form a square of side 3
    for i in range(rows - 2):
        for j in range(cols - 2):
            # Check for squares of increasing size starting from 3x3
            for size in range(3, min(rows - i, cols - j) + 1):
                # Check if all elements in the square are the same
                char = matrix[i][j]
                is_square = True
                for x in range(i, i + size):
                    for y in range(j, j + size):
                        if matrix[x][y] != char:
                            is_square = False
                            break
                    if not is_square:
                        break
                if is_square:
                    # Return the bottom-right corner of the square
                    return [i + size - 1, j + size - 1]

# Define the matrix
matrix = [
    ['J', 'K', 'S', 'P', 'X', 'L'],
    ['Q', 'V', 'G', 'W', 'W', 'W'],
    ['T', 'U', 'S', 'W', 'W', 'W'],
    ['N', 'X', 'N', 'W', 'W', 'W'],
    ['E', 'I', 'L', 'H', 'A', 'B'],
    ['R', 'Z', 'G', 'U', 'D', 'V'],
    ['J', 'S', 'B', 'F', 'U', 'X'],
    ['Y', 'K', 'O', 'V', 'U', 'E'],
    ['Z', 'F', 'X', 'R', 'I', 'N']
]

# Find and print the bottom-right corner of the square
print(find_square(matrix))
def find_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    
    # We need at least a 3x3 matrix to form a square of side 3
    for r in range(rows - 2):
        for c in range(cols - 2):
            # Check for squares of increasing size starting from 3x3
            max_size = min(rows - r, cols - c)
            for size in range(3, max_size + 1):
                # Check if all elements in the square are the same
                char = matrix[r][c]
                is_square = True
                for i in range(size):
                    for j in range(size):
                        if matrix[r + i][c + j] != char:
                            is_square = False
                            break
                    if not is_square:
                        break
                if is_square:
                    # Return the bottom-right corner of the square
                    return [r + size - 1, c + size - 1]

# Define the matrix
matrix = [
    ['C', 'Z', 'Z', 'C', 'Z', 'U', 'K', 'V', 'M', 'D', 'D', 'F', 'O', 'J'],
    ['W', 'V', 'F', 'S', 'A', 'U', 'G', 'J', 'X', 'P', 'D', 'R', 'J', 'D'],
    ['H', 'H', 'T', 'G', 'Q', 'X', 'F', 'D', 'Z', 'T', 'T', 'T', 'T', 'D'],
    ['J', 'I', 'M', 'Y', 'B', 'T', 'W', 'X', 'K', 'T', 'T', 'T', 'T', 'H'],
    ['D', 'U', 'H', 'E', 'F', 'P', 'L', 'I', 'I', 'T', 'T', 'T', 'T', 'O'],
    ['S', 'A', 'Y', 'W', 'M', 'V', 'M', 'M', 'T', 'T', 'T', 'T', 'T', 'N'],
    ['T', 'K', 'B', 'K', 'M', 'J', 'Z', 'V', 'T', 'S', 'A', 'T', 'E', 'M'],
    ['G', 'X', 'N', 'E', 'Q', 'K', 'W', 'R', 'B', 'E', 'I', 'W', 'Y', 'J'],
    ['Z', 'T', 'I', 'I', 'M', 'M', 'R', 'E', 'A', 'R', 'Q', 'W', 'F', 'V'],
    ['K', 'B', 'K', 'D', 'J', 'F', 'C', 'C', 'X', 'Q', 'P', 'M', 'E', 'O'],
    ['H', 'D', 'Z', 'K', 'T', 'Q', 'K', 'S', 'T', 'L', 'W', 'H', 'D', 'Z']
]

# Find the bottom-right corner of the square
result = find_square(matrix)
print(result)
def find_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    
    # We need at least a 3x3 area to form a square of side 3
    for i in range(rows - 2):
        for j in range(cols - 2):
            # Check for squares of increasing size starting from 3
            max_size = min(rows - i, cols - j)
            for size in range(3, max_size + 1):
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
    ['Y', 'V', 'F', 'S', 'F', 'A', 'T', 'P'],
    ['M', 'X', 'Z', 'P', 'Z', 'D', 'O', 'T'],
    ['U', 'O', 'O', 'X', 'U', 'M', 'E', 'Y'],
    ['E', 'T', 'A', 'G', 'Z', 'Q', 'A', 'L'],
    ['F', 'R', 'R', 'H', 'R', 'P', 'R', 'P'],
    ['G', 'T', 'K', 'T', 'Z', 'W', 'W', 'W'],
    ['M', 'S', 'E', 'L', 'Z', 'W', 'W', 'W'],
    ['E', 'G', 'Y', 'W', 'S', 'W', 'W', 'W'],
    ['Z', 'H', 'T', 'K', 'F', 'F', 'Y', 'S'],
    ['O', 'R', 'V', 'H', 'V', 'S', 'S', 'V'],
    ['F', 'J', 'X', 'J', 'X', 'X', 'B', 'R'],
    ['C', 'O', 'G', 'T', 'J', 'G', 'L', 'K'],
    ['A', 'A', 'R', 'R', 'D', 'L', 'R', 'I']
]

# Find and print the bottom-right corner of the square
print(find_square(matrix))
def find_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    
    for i in range(rows):
        for j in range(cols):
            # Check for squares of size 3x3 or larger
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
    ['G', 'R', 'Z', 'K', 'I', 'V'],
    ['C', 'J', 'H', 'B', 'Y', 'K'],
    ['B', 'L', 'L', 'Q', 'U', 'X'],
    ['U', 'O', 'R', 'Q', 'T', 'R'],
    ['Z', 'E', 'L', 'O', 'O', 'O'],
    ['X', 'J', 'J', 'O', 'O', 'O'],
    ['C', 'Y', 'M', 'O', 'O', 'O']
]

# Find and print the bottom-right corner of the square
print(find_square(matrix))
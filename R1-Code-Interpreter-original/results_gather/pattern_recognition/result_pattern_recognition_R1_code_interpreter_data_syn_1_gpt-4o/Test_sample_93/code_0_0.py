def find_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    
    # Start checking for squares of size 3x3 and larger
    for size in range(3, min(rows, cols) + 1):
        for i in range(rows - size + 1):
            for j in range(cols - size + 1):
                # Check if all characters in the square are the same
                char = matrix[i][j]
                is_square = True
                for x in range(size):
                    for y in range(size):
                        if matrix[i + x][j + y] != char:
                            is_square = False
                            break
                    if not is_square:
                        break
                if is_square:
                    # Return the bottom-right corner of the square
                    return [i + size - 1, j + size - 1]

# Define the matrix
matrix = [
    ['A', 'S', 'U', 'A', 'X', 'S', 'D', 'F', 'L'],
    ['W', 'X', 'M', 'U', 'N', 'C', 'X', 'F', 'E'],
    ['E', 'B', 'B', 'P', 'S', 'N', 'Y', 'Z', 'N'],
    ['F', 'T', 'J', 'O', 'V', 'O', 'D', 'P', 'E'],
    ['F', 'Y', 'Q', 'I', 'R', 'P', 'N', 'L', 'Y'],
    ['S', 'W', 'J', 'Y', 'X', 'C', 'P', 'V', 'P'],
    ['K', 'O', 'Q', 'D', 'C', 'G', 'Z', 'B', 'U'],
    ['F', 'Y', 'L', 'X', 'F', 'E', 'Z', 'S', 'H'],
    ['E', 'W', 'O', 'J', 'R', 'X', 'O', 'E', 'X'],
    ['O', 'H', 'H', 'H', 'S', 'U', 'V', 'X', 'G'],
    ['Y', 'H', 'H', 'H', 'E', 'S', 'X', 'I', 'J'],
    ['H', 'H', 'H', 'H', 'I', 'I', 'R', 'W', 'U'],
    ['E', 'O', 'B', 'M', 'B', 'Y', 'F', 'T', 'K'],
    ['B', 'H', 'Z', 'I', 'Z', 'O', 'X', 'A', 'V']
]

# Find and print the bottom-right corner of the square
print(find_square(matrix))
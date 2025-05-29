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
    ['T', 'N', 'G', 'F', 'I', 'M', 'O'],
    ['I', 'U', 'G', 'I', 'C', 'S', 'V'],
    ['I', 'N', 'U', 'Z', 'R', 'Y', 'C'],
    ['V', 'V', 'Y', 'A', 'G', 'O', 'R'],
    ['Y', 'N', 'Z', 'Y', 'D', 'M', 'V'],
    ['P', 'W', 'S', 'X', 'C', 'E', 'K'],
    ['Q', 'B', 'W', 'P', 'O', 'A', 'D'],
    ['N', 'Z', 'C', 'P', 'T', 'C', 'N'],
    ['J', 'A', 'M', 'Z', 'H', 'E', 'J'],
    ['S', 'Q', 'O', 'Z', 'B', 'P', 'B'],
    ['X', 'X', 'X', 'U', 'C', 'C', 'O'],
    ['X', 'X', 'X', 'O', 'J', 'U', 'R'],
    ['X', 'X', 'X', 'O', 'L', 'R', 'H']
]

# Find and print the bottom-right corner of the square
print(find_square(matrix))
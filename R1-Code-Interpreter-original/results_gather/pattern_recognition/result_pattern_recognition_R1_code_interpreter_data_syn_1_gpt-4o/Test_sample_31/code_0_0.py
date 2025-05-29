def find_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    
    # Iterate over each possible top-left corner of a square
    for i in range(rows - 2):
        for j in range(cols - 2):
            # Check for squares of increasing size starting from 3
            for size in range(3, min(rows - i, cols - j) + 1):
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
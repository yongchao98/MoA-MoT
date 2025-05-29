def find_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    
    # Iterate over each possible top-left corner of a square
    for i in range(rows - 2):
        for j in range(cols - 2):
            # Check for squares of increasing size
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
    ['M', 'Y', 'E', 'O', 'Q', 'B', 'O', 'N'],
    ['P', 'P', 'D', 'M', 'B', 'F', 'X', 'P'],
    ['J', 'E', 'X', 'O', 'O', 'O', 'O', 'O'],
    ['V', 'G', 'C', 'O', 'O', 'O', 'O', 'W'],
    ['H', 'Z', 'K', 'O', 'O', 'O', 'O', 'O'],
    ['J', 'V', 'L', 'O', 'O', 'O', 'O', 'E'],
    ['K', 'Y', 'N', 'Z', 'Q', 'M', 'F', 'B'],
    ['F', 'R', 'T', 'N', 'W', 'X', 'O', 'T'],
    ['T', 'A', 'D', 'X', 'U', 'N', 'U', 'U'],
    ['X', 'B', 'X', 'E', 'Y', 'M', 'R', 'W']
]

# Find and print the bottom-right corner of the square
print(find_square(matrix))
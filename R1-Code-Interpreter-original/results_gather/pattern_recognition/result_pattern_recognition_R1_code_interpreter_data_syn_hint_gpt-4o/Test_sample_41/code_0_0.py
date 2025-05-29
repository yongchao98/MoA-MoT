def find_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    
    for i in range(rows):
        for j in range(cols):
            # Check for squares of size 3x3 or larger
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
    ['T', 'E', 'E', 'E', 'X', 'V', 'M', 'P', 'Z', 'I', 'V', 'L'],
    ['P', 'E', 'E', 'E', 'F', 'H', 'E', 'M', 'X', 'Q', 'P', 'W'],
    ['R', 'E', 'E', 'E', 'C', 'O', 'Q', 'S', 'O', 'T', 'D', 'S'],
    ['N', 'I', 'B', 'O', 'X', 'T', 'D', 'F', 'I', 'V', 'A', 'T'],
    ['U', 'T', 'X', 'P', 'E', 'Y', 'U', 'Y', 'P', 'F', 'N', 'F'],
    ['H', 'B', 'Z', 'I', 'M', 'N', 'X', 'Y', 'A', 'G', 'J', 'W']
]

# Find and print the bottom-right corner of the square
print(find_square(matrix))
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
    ['Q', 'L', 'H', 'G', 'V', 'I', 'O', 'V', 'H', 'N', 'O', 'O', 'X', 'K', 'I'],
    ['N', 'W', 'H', 'X', 'B', 'E', 'K', 'C', 'W', 'T', 'A', 'J', 'I', 'O', 'D'],
    ['X', 'Q', 'B', 'V', 'S', 'X', 'X', 'I', 'D', 'I', 'Y', 'C', 'T', 'N', 'L'],
    ['O', 'U', 'I', 'N', 'Y', 'W', 'C', 'I', 'I', 'I', 'I', 'P', 'P', 'Q', 'H'],
    ['M', 'R', 'N', 'U', 'X', 'O', 'Y', 'U', 'I', 'I', 'I', 'J', 'J', 'Y', 'T'],
    ['A', 'T', 'E', 'U', 'W', 'Y', 'K', 'Z', 'I', 'I', 'I', 'K', 'L', 'W', 'S']
]

# Find and print the bottom-right corner of the square
print(find_square(matrix))
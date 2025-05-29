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
    ['X', 'R', 'C', 'H', 'H', 'E', 'Y', 'Y', 'Y', 'E', 'N', 'A'],
    ['D', 'L', 'W', 'I', 'Q', 'R', 'Y', 'Y', 'Y', 'V', 'G', 'P'],
    ['Y', 'D', 'R', 'J', 'L', 'D', 'Y', 'Y', 'Y', 'R', 'Q', 'B'],
    ['A', 'D', 'C', 'S', 'J', 'W', 'S', 'G', 'O', 'T', 'Q', 'V'],
    ['G', 'Y', 'J', 'O', 'G', 'U', 'Z', 'C', 'H', 'E', 'H', 'G'],
    ['I', 'H', 'X', 'N', 'R', 'M', 'Q', 'Z', 'G', 'J', 'V', 'U'],
    ['O', 'H', 'L', 'O', 'Y', 'A', 'W', 'I', 'M', 'W', 'F', 'J']
]

# Find and print the bottom-right corner of the square
print(find_square(matrix))
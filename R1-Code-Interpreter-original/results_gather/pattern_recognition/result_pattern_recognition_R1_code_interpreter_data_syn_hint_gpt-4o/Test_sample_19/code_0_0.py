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
    ['Q', 'L', 'H', 'E', 'X', 'Q', 'Y', 'D', 'Y', 'B', 'E', 'F'],
    ['L', 'I', 'A', 'M', 'B', 'G', 'X', 'C', 'N', 'J', 'V', 'W'],
    ['C', 'Q', 'W', 'N', 'F', 'A', 'B', 'U', 'O', 'P', 'B', 'X'],
    ['F', 'O', 'K', 'I', 'L', 'A', 'Z', 'F', 'L', 'L', 'L', 'G'],
    ['I', 'V', 'S', 'Q', 'N', 'T', 'M', 'C', 'L', 'L', 'L', 'E'],
    ['G', 'T', 'V', 'P', 'K', 'H', 'G', 'L', 'L', 'L', 'L', 'I']
]

# Find and print the bottom-right corner of the square
print(find_square(matrix))
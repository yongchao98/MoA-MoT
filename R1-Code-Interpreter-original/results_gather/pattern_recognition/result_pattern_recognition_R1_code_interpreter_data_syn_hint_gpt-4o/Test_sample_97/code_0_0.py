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
    ['Y', 'L', 'I', 'F', 'T', 'Q', 'D', 'Y', 'F', 'I', 'V', 'X', 'L', 'V'],
    ['F', 'S', 'H', 'W', 'Y', 'J', 'X', 'R', 'I', 'T', 'I', 'R', 'B', 'H'],
    ['T', 'T', 'X', 'X', 'R', 'P', 'G', 'N', 'I', 'F', 'Q', 'T', 'N', 'A'],
    ['D', 'T', 'C', 'N', 'U', 'V', 'G', 'E', 'D', 'E', 'J', 'Z', 'M', 'R'],
    ['K', 'S', 'B', 'R', 'S', 'N', 'V', 'W', 'S', 'D', 'K', 'Y', 'W', 'J'],
    ['M', 'H', 'P', 'E', 'A', 'D', 'M', 'V', 'X', 'N', 'E', 'E', 'E', 'P'],
    ['D', 'Y', 'K', 'V', 'Q', 'H', 'N', 'U', 'A', 'X', 'E', 'E', 'E', 'A'],
    ['Z', 'G', 'D', 'R', 'D', 'M', 'F', 'Z', 'M', 'B', 'E', 'E', 'E', 'W']
]

# Find and print the bottom-right corner of the square
print(find_square(matrix))
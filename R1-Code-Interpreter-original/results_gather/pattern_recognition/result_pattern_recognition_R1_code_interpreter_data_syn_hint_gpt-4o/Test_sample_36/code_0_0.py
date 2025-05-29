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
    ['S', 'P', 'W', 'I', 'W', 'H', 'N'],
    ['U', 'L', 'S', 'P', 'N', 'N', 'N'],
    ['M', 'V', 'H', 'H', 'N', 'N', 'N'],
    ['I', 'R', 'I', 'F', 'N', 'N', 'N'],
    ['W', 'M', 'I', 'D', 'F', 'P', 'S'],
    ['T', 'H', 'X', 'M', 'U', 'K', 'C'],
    ['Y', 'A', 'W', 'N', 'L', 'P', 'H'],
    ['O', 'V', 'I', 'S', 'N', 'Q', 'T'],
    ['T', 'K', 'O', 'S', 'E', 'L', 'I'],
    ['B', 'Y', 'Y', 'G', 'N', 'B', 'S'],
    ['U', 'Y', 'E', 'E', 'M', 'W', 'J'],
    ['W', 'O', 'C', 'C', 'V', 'L', 'W']
]

# Find and print the bottom-right corner of the square
print(find_square(matrix))
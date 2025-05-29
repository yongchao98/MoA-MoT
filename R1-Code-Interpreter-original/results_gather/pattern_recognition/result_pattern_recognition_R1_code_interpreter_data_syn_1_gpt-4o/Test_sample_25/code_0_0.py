def find_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    
    # We need at least a 3x3 matrix to form a square of size 3
    for size in range(3, min(rows, cols) + 1):
        for i in range(rows - size + 1):
            for j in range(cols - size + 1):
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
    ['N', 'H', 'M', 'K', 'H', 'L', 'X', 'G', 'S', 'R'],
    ['A', 'L', 'Z', 'D', 'Z', 'S', 'S', 'R', 'G', 'X'],
    ['N', 'O', 'P', 'I', 'I', 'I', 'I', 'B', 'I', 'H'],
    ['R', 'R', 'Y', 'I', 'I', 'I', 'I', 'H', 'D', 'V'],
    ['J', 'I', 'X', 'I', 'I', 'I', 'I', 'D', 'Y', 'I'],
    ['T', 'N', 'G', 'I', 'I', 'I', 'I', 'Z', 'L', 'J'],
    ['C', 'V', 'E', 'T', 'D', 'P', 'Q', 'T', 'Y', 'B'],
    ['C', 'J', 'A', 'N', 'S', 'S', 'J', 'S', 'O', 'E']
]

# Find and print the bottom-right corner of the square
print(find_square(matrix))
def find_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    
    # We need to find the largest square with side length at least 3
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
    ['Z', 'X', 'O', 'F', 'I', 'Y', 'C', 'K', 'K', 'A', 'I'],
    ['B', 'G', 'Z', 'H', 'X', 'X', 'J', 'A', 'P', 'V', 'B'],
    ['K', 'T', 'B', 'A', 'F', 'W', 'T', 'N', 'Y', 'X', 'T'],
    ['H', 'F', 'L', 'Q', 'W', 'J', 'S', 'I', 'Q', 'H', 'L'],
    ['Y', 'I', 'Q', 'R', 'Z', 'F', 'K', 'U', 'L', 'R', 'C'],
    ['X', 'F', 'K', 'J', 'J', 'J', 'J', 'N', 'X', 'H', 'K'],
    ['D', 'L', 'U', 'J', 'J', 'J', 'J', 'X', 'L', 'Z', 'E'],
    ['S', 'H', 'G', 'J', 'J', 'J', 'J', 'V', 'Z', 'R', 'Y'],
    ['O', 'M', 'X', 'J', 'J', 'J', 'J', 'W', 'M', 'I', 'R'],
    ['P', 'O', 'C', 'O', 'I', 'N', 'F', 'S', 'M', 'V', 'Y']
]

# Find and print the bottom-right corner of the square
print(find_square(matrix))
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
    ['I', 'P', 'U', 'S', 'N', 'S', 'A'],
    ['L', 'C', 'S', 'P', 'E', 'X', 'M'],
    ['I', 'F', 'Z', 'C', 'A', 'K', 'S'],
    ['Y', 'E', 'C', 'H', 'B', 'M', 'Y'],
    ['R', 'C', 'F', 'C', 'C', 'Q', 'F'],
    ['R', 'J', 'F', 'M', 'X', 'X', 'K'],
    ['I', 'I', 'I', 'L', 'U', 'Q', 'M'],
    ['I', 'I', 'I', 'I', 'P', 'K', 'F'],
    ['I', 'I', 'I', 'M', 'U', 'P', 'D'],
    ['Q', 'L', 'F', 'A', 'R', 'W', 'K'],
    ['X', 'P', 'A', 'X', 'I', 'U', 'Z'],
    ['T', 'Z', 'U', 'X', 'P', 'B', 'F']
]

# Find and print the bottom-right corner of the square
print(find_square(matrix))
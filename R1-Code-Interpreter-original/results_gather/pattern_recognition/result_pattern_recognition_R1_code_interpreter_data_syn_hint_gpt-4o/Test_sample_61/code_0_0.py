def find_largest_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    
    # Start with the largest possible square size and work downwards
    for size in range(min(rows, cols), 2, -1):
        for r in range(rows - size + 1):
            for c in range(cols - size + 1):
                # Check if all characters in the square are the same
                char = matrix[r][c]
                is_square = True
                for i in range(size):
                    for j in range(size):
                        if matrix[r + i][c + j] != char:
                            is_square = False
                            break
                    if not is_square:
                        break
                if is_square:
                    # Return the bottom-right corner of the square
                    return [r + size - 1, c + size - 1]

# Define the matrix
matrix = [
    ['Z', 'W', 'P', 'T', 'I', 'V', 'J', 'I', 'Y'],
    ['C', 'C', 'J', 'T', 'M', 'B', 'I', 'F', 'Y'],
    ['Q', 'A', 'F', 'F', 'F', 'F', 'F', 'M', 'Q'],
    ['J', 'L', 'C', 'F', 'F', 'F', 'F', 'J', 'B'],
    ['U', 'V', 'M', 'F', 'F', 'F', 'F', 'D', 'N'],
    ['M', 'G', 'I', 'F', 'F', 'F', 'F', 'S', 'Q'],
    ['G', 'E', 'T', 'C', 'N', 'M', 'S', 'G', 'U'],
    ['R', 'S', 'U', 'X', 'D', 'D', 'O', 'P', 'A'],
    ['H', 'B', 'I', 'R', 'S', 'K', 'B', 'K', 'Z'],
    ['U', 'N', 'W', 'Z', 'I', 'M', 'J', 'Z', 'Q'],
    ['K', 'I', 'D', 'J', 'H', 'M', 'C', 'Y', 'R']
]

# Find and print the bottom-right corner of the largest square
print(find_largest_square(matrix))
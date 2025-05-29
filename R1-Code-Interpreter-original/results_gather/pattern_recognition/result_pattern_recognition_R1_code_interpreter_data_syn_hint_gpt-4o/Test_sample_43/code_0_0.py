def find_largest_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    max_size = min(rows, cols)
    
    # Start checking from the largest possible square size
    for size in range(max_size, 2, -1):
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
                    return [i + size - 1, j + size - 1]

# Define the matrix
matrix = [
    ['O', 'D', 'U', 'U', 'B', 'X', 'L', 'D', 'W', 'X', 'N', 'N', 'U', 'T', 'V'],
    ['U', 'D', 'J', 'V', 'F', 'T', 'W', 'F', 'L', 'D', 'O', 'S', 'T', 'Z', 'U'],
    ['A', 'W', 'B', 'B', 'Z', 'M', 'I', 'J', 'F', 'W', 'V', 'M', 'Z', 'P', 'V'],
    ['T', 'L', 'D', 'N', 'P', 'E', 'Y', 'C', 'V', 'U', 'E', 'X', 'V', 'V', 'Z'],
    ['E', 'E', 'E', 'E', 'E', 'Z', 'W', 'A', 'A', 'Q', 'G', 'U', 'B', 'K', 'X'],
    ['E', 'E', 'E', 'E', 'E', 'U', 'G', 'J', 'U', 'G', 'N', 'H', 'Z', 'Q', 'H'],
    ['E', 'E', 'E', 'E', 'E', 'Q', 'X', 'S', 'W', 'Z', 'F', 'F', 'P', 'O', 'C'],
    ['E', 'E', 'E', 'E', 'E', 'V', 'Q', 'Z', 'D', 'W', 'X', 'L', 'Q', 'C', 'H'],
    ['E', 'E', 'E', 'E', 'E', 'E', 'A', 'L', 'Z', 'M', 'Y', 'E', 'D', 'N', 'B'],
    ['V', 'K', 'F', 'Y', 'R', 'X', 'K', 'G', 'A', 'V', 'X', 'M', 'Y', 'T', 'R'],
    ['S', 'J', 'J', 'N', 'P', 'O', 'W', 'R', 'J', 'O', 'P', 'R', 'K', 'M', 'R'],
    ['H', 'P', 'D', 'J', 'X', 'C', 'T', 'V', 'M', 'W', 'C', 'I', 'W', 'T', 'M']
]

# Find and print the bottom-right corner of the largest square
print(find_largest_square(matrix))
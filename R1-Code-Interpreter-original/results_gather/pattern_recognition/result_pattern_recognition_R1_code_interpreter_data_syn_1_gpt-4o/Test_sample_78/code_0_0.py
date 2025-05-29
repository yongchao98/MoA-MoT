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
    ['V', 'J', 'Q', 'W', 'S', 'X', 'C', 'B'],
    ['T', 'X', 'R', 'K', 'K', 'Z', 'N', 'E'],
    ['J', 'T', 'X', 'R', 'A', 'S', 'L', 'X'],
    ['J', 'H', 'H', 'J', 'M', 'P', 'P', 'P'],
    ['Q', 'T', 'I', 'F', 'H', 'P', 'P', 'P'],
    ['A', 'K', 'T', 'X', 'R', 'P', 'P', 'P'],
    ['P', 'N', 'J', 'N', 'T', 'G', 'K', 'O']
]

# Find and print the bottom-right corner of the square
print(find_square(matrix))
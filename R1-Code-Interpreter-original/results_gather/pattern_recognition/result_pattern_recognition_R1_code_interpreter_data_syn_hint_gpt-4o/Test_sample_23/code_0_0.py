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
    ['N', 'Q', 'S', 'F', 'Z', 'W', 'T', 'N', 'M', 'S', 'H', 'G'],
    ['V', 'R', 'B', 'K', 'X', 'F', 'D', 'N', 'X', 'F', 'Q', 'S'],
    ['Z', 'B', 'N', 'U', 'L', 'I', 'V', 'C', 'V', 'G', 'W', 'B'],
    ['B', 'H', 'R', 'X', 'J', 'O', 'U', 'Q', 'E', 'O', 'T', 'G'],
    ['M', 'R', 'S', 'G', 'C', 'P', 'B', 'A', 'G', 'X', 'H', 'M'],
    ['A', 'R', 'K', 'N', 'G', 'G', 'P', 'A', 'E', 'R', 'Z', 'J'],
    ['J', 'A', 'K', 'Y', 'A', 'Z', 'I', 'V', 'E', 'R', 'B', 'X'],
    ['B', 'L', 'H', 'P', 'D', 'V', 'Y', 'X', 'N', 'B', 'G', 'K'],
    ['Q', 'U', 'B', 'G', 'S', 'Q', 'W', 'M', 'V', 'U', 'R', 'V'],
    ['B', 'C', 'L', 'L', 'H', 'X', 'K', 'B', 'U', 'H', 'T', 'N'],
    ['N', 'P', 'N', 'E', 'J', 'N', 'X', 'Z', 'J', 'S', 'L', 'K'],
    ['T', 'Q', 'J', 'D', 'T', 'W', 'N', 'K', 'Q', 'F', 'E', 'X'],
    ['K', 'X', 'I', 'C', 'O', 'L', 'S', 'M', 'M', 'M', 'B', 'Y'],
    ['X', 'J', 'W', 'Q', 'K', 'L', 'G', 'M', 'M', 'M', 'C', 'M'],
    ['V', 'D', 'T', 'N', 'X', 'L', 'I', 'M', 'M', 'M', 'T', 'A']
]

# Find and print the bottom-right corner of the square
print(find_square(matrix))
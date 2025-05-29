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
    ['E', 'Z', 'I', 'I', 'I', 'M', 'Y', 'V', 'F', 'O', 'X', 'U', 'Q', 'Y'],
    ['B', 'Y', 'I', 'I', 'I', 'Z', 'Y', 'F', 'X', 'G', 'G', 'Z', 'M', 'I'],
    ['L', 'C', 'I', 'I', 'I', 'J', 'L', 'M', 'S', 'U', 'G', 'E', 'T', 'X'],
    ['B', 'S', 'T', 'R', 'E', 'U', 'S', 'A', 'O', 'K', 'M', 'B', 'X', 'A'],
    ['N', 'F', 'T', 'S', 'C', 'A', 'E', 'C', 'W', 'Q', 'V', 'T', 'A', 'H'],
    ['T', 'D', 'V', 'U', 'H', 'M', 'E', 'H', 'S', 'Q', 'Q', 'R', 'Z', 'O']
]

# Find and print the bottom-right corner of the square
print(find_square(matrix))
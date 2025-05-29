def find_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    
    # Iterate over possible top-left corners of squares
    for i in range(rows - 2):
        for j in range(cols - 2):
            # Check for squares of increasing size starting from 3x3
            for size in range(3, min(rows - i, cols - j) + 1):
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
    ['O', 'G', 'U', 'F', 'B', 'F', 'C'],
    ['M', 'O', 'M', 'K', 'M', 'F', 'D'],
    ['X', 'V', 'N', 'C', 'I', 'P', 'M'],
    ['U', 'Z', 'U', 'G', 'P', 'U', 'A'],
    ['E', 'C', 'Q', 'P', 'C', 'C', 'C'],
    ['G', 'P', 'A', 'E', 'C', 'C', 'C'],
    ['U', 'D', 'Y', 'D', 'C', 'C', 'C'],
    ['R', 'S', 'G', 'W', 'F', 'A', 'Z']
]

# Find and print the bottom-right corner of the square
print(find_square(matrix))
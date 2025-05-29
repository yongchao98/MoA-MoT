def find_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    
    # Iterate over each possible top-left corner of a square
    for i in range(rows - 2):
        for j in range(cols - 2):
            # Check for squares of increasing size starting from 3x3
            max_size = min(rows - i, cols - j)
            for size in range(3, max_size + 1):
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
                    # Record the bottom-right corner of the square
                    bottom_right = (i + size - 1, j + size - 1)
                    return bottom_right

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

# Find the bottom-right corner of the square
result = find_square(matrix)
print(result)
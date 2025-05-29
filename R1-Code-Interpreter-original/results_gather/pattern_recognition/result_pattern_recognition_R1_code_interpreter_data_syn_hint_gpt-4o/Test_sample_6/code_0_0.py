def find_largest_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    max_size = 0
    bottom_right_corner = None

    # Iterate over each possible top-left corner of a square
    for i in range(rows):
        for j in range(cols):
            # Check for squares of increasing size
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
                
                # If a valid square is found, update the largest square info
                if is_square and size > max_size:
                    max_size = size
                    bottom_right_corner = (i + size - 1, j + size - 1)

    return bottom_right_corner

# Define the matrix
matrix = [
    ['L', 'W', 'F', 'E', 'O', 'J'],
    ['V', 'J', 'W', 'O', 'V', 'E'],
    ['F', 'Y', 'F', 'W', 'V', 'T'],
    ['Z', 'Y', 'Z', 'U', 'Q', 'K'],
    ['Y', 'Y', 'K', 'W', 'F', 'G'],
    ['J', 'O', 'S', 'A', 'H', 'L'],
    ['B', 'G', 'E', 'I', 'K', 'S'],
    ['C', 'S', 'A', 'A', 'A', 'U'],
    ['Y', 'A', 'A', 'A', 'A', 'N'],
    ['U', 'Y', 'A', 'A', 'A', 'G']
]

# Find the bottom-right corner of the largest square
result = find_largest_square(matrix)
print(result)
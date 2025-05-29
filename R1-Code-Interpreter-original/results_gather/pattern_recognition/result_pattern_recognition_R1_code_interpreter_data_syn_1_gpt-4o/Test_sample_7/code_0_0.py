def find_largest_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    max_size = 0
    bottom_right_corner = (-1, -1)

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
                
                # If a valid square is found, update the largest square
                if is_square and size > max_size:
                    max_size = size
                    bottom_right_corner = (i + size - 1, j + size - 1)

    return bottom_right_corner

# Define the matrix
matrix = [
    ['C', 'S', 'A', 'F', 'O', 'U', 'I', 'T', 'C'],
    ['Y', 'Z', 'N', 'H', 'Y', 'B', 'B', 'G', 'O'],
    ['S', 'J', 'C', 'X', 'H', 'C', 'T', 'Q', 'I'],
    ['K', 'S', 'W', 'F', 'C', 'S', 'P', 'G', 'X'],
    ['G', 'G', 'G', 'G', 'Y', 'Q', 'A', 'Y', 'G'],
    ['G', 'G', 'G', 'G', 'M', 'S', 'M', 'K', 'M'],
    ['G', 'G', 'G', 'G', 'S', 'O', 'P', 'G', 'H'],
    ['G', 'G', 'G', 'G', 'G', 'A', 'Y', 'X', 'X'],
    ['D', 'X', 'V', 'U', 'X', 'R', 'A', 'L', 'T'],
    ['T', 'B', 'M', 'F', 'J', 'B', 'H', 'W', 'Q'],
    ['B', 'N', 'Z', 'K', 'Y', 'T', 'D', 'E', 'Q'],
    ['Y', 'E', 'T', 'Z', 'T', 'M', 'E', 'Z', 'K'],
    ['W', 'A', 'R', 'T', 'V', 'X', 'Z', 'H', 'X']
]

# Find the bottom-right corner of the largest square
result = find_largest_square(matrix)
print(result)
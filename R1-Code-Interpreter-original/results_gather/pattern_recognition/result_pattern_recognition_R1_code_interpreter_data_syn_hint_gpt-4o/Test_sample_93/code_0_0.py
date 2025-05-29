def find_largest_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    max_size = 0
    bottom_right_corner = (-1, -1)

    # Iterate over each possible top-left corner of a square
    for i in range(rows):
        for j in range(cols):
            # Try to form squares of increasing size
            for size in range(3, min(rows - i, cols - j) + 1):
                # Check if the square of this size is valid
                char = matrix[i][j]
                valid_square = True
                # Check the bottom and right boundaries of the square
                for k in range(size):
                    if matrix[i + size - 1][j + k] != char or matrix[i + k][j + size - 1] != char:
                        valid_square = False
                        break
                if valid_square:
                    if size > max_size:
                        max_size = size
                        bottom_right_corner = (i + size - 1, j + size - 1)

    return bottom_right_corner

matrix = [
    ['A', 'S', 'U', 'A', 'X', 'S', 'D', 'F', 'L'],
    ['W', 'X', 'M', 'U', 'N', 'C', 'X', 'F', 'E'],
    ['E', 'B', 'B', 'P', 'S', 'N', 'Y', 'Z', 'N'],
    ['F', 'T', 'J', 'O', 'V', 'O', 'D', 'P', 'E'],
    ['F', 'Y', 'Q', 'I', 'R', 'P', 'N', 'L', 'Y'],
    ['S', 'W', 'J', 'Y', 'X', 'C', 'P', 'V', 'P'],
    ['K', 'O', 'Q', 'D', 'C', 'G', 'Z', 'B', 'U'],
    ['F', 'Y', 'L', 'X', 'F', 'E', 'Z', 'S', 'H'],
    ['E', 'W', 'O', 'J', 'R', 'X', 'O', 'E', 'X'],
    ['O', 'H', 'H', 'H', 'S', 'U', 'V', 'X', 'G'],
    ['Y', 'H', 'H', 'H', 'E', 'S', 'X', 'I', 'J'],
    ['H', 'H', 'H', 'H', 'I', 'I', 'R', 'W', 'U'],
    ['E', 'O', 'B', 'M', 'B', 'Y', 'F', 'T', 'K'],
    ['B', 'H', 'Z', 'I', 'Z', 'O', 'X', 'A', 'V']
]

result = find_largest_square(matrix)
print(result)
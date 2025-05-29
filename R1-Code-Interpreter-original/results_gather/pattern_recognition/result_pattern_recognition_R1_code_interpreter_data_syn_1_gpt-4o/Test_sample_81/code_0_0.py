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
                # Check the right and bottom boundaries of the square
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
    ['Q', 'L', 'H', 'G', 'V', 'I', 'O', 'V', 'H', 'N', 'O', 'O', 'X', 'K', 'I'],
    ['N', 'W', 'H', 'X', 'B', 'E', 'K', 'C', 'W', 'T', 'A', 'J', 'I', 'O', 'D'],
    ['X', 'Q', 'B', 'V', 'S', 'X', 'X', 'I', 'D', 'I', 'Y', 'C', 'T', 'N', 'L'],
    ['O', 'U', 'I', 'N', 'Y', 'W', 'C', 'I', 'I', 'I', 'I', 'P', 'P', 'Q', 'H'],
    ['M', 'R', 'N', 'U', 'X', 'O', 'Y', 'U', 'I', 'I', 'I', 'J', 'J', 'Y', 'T'],
    ['A', 'T', 'E', 'U', 'W', 'Y', 'K', 'Z', 'I', 'I', 'I', 'K', 'L', 'W', 'S']
]

result = find_largest_square(matrix)
print(result)
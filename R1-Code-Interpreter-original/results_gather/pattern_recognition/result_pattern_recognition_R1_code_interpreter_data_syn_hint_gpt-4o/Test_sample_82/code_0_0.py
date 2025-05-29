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
                # Check if a square of this size can be formed
                char = matrix[i][j]
                valid_square = True
                # Check the right and bottom boundaries
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
    ['P', 'S', 'K', 'Y', 'A', 'O', 'X', 'B', 'H', 'C'],
    ['T', 'Y', 'Z', 'V', 'T', 'M', 'K', 'X', 'A', 'T'],
    ['D', 'P', 'R', 'O', 'Z', 'K', 'K', 'K', 'Y', 'J'],
    ['R', 'X', 'K', 'W', 'M', 'K', 'K', 'K', 'N', 'G'],
    ['B', 'V', 'B', 'T', 'D', 'K', 'K', 'K', 'K', 'I'],
    ['B', 'M', 'I', 'D', 'X', 'M', 'P', 'E', 'A', 'K']
]

result = find_largest_square(matrix)
print(result)
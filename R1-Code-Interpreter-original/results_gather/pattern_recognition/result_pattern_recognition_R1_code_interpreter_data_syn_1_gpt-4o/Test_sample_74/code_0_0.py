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
                    if matrix[i + k][j + size - 1] != char or matrix[i + size - 1][j + k] != char:
                        valid_square = False
                        break
                if valid_square:
                    if size > max_size:
                        max_size = size
                        bottom_right_corner = (i + size - 1, j + size - 1)

    return bottom_right_corner

matrix = [
    ['J', 'A', 'L', 'T', 'O', 'H'],
    ['Z', 'U', 'B', 'V', 'M', 'Z'],
    ['P', 'W', 'U', 'V', 'M', 'F'],
    ['S', 'R', 'Q', 'C', 'P', 'S'],
    ['N', 'Y', 'O', 'Q', 'T', 'M'],
    ['H', 'A', 'W', 'W', 'W', 'T'],
    ['K', 'S', 'W', 'W', 'W', 'P'],
    ['F', 'O', 'W', 'W', 'W', 'V'],
    ['L', 'Q', 'X', 'W', 'M', 'S'],
    ['U', 'O', 'O', 'W', 'Q', 'M'],
    ['I', 'I', 'Q', 'V', 'Y', 'M'],
    ['R', 'S', 'R', 'U', 'A', 'H']
]

result = find_largest_square(matrix)
print(result)
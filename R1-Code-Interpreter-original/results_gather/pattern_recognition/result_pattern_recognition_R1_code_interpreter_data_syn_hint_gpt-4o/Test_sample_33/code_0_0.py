def find_largest_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    max_size = 0
    bottom_right_corner = (-1, -1)

    # Iterate over each possible starting point for a square
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
    ['C', 'Z', 'Z', 'C', 'Z', 'U', 'K', 'V', 'M', 'D', 'D', 'F', 'O', 'J'],
    ['W', 'V', 'F', 'S', 'A', 'U', 'G', 'J', 'X', 'P', 'D', 'R', 'J', 'D'],
    ['H', 'H', 'T', 'G', 'Q', 'X', 'F', 'D', 'Z', 'T', 'T', 'T', 'T', 'D'],
    ['J', 'I', 'M', 'Y', 'B', 'T', 'W', 'X', 'K', 'T', 'T', 'T', 'T', 'H'],
    ['D', 'U', 'H', 'E', 'F', 'P', 'L', 'I', 'I', 'T', 'T', 'T', 'T', 'O'],
    ['S', 'A', 'Y', 'W', 'M', 'V', 'M', 'M', 'T', 'T', 'T', 'T', 'T', 'N'],
    ['T', 'K', 'B', 'K', 'M', 'J', 'Z', 'V', 'T', 'S', 'A', 'T', 'E', 'M'],
    ['G', 'X', 'N', 'E', 'Q', 'K', 'W', 'R', 'B', 'E', 'I', 'W', 'Y', 'J'],
    ['Z', 'T', 'I', 'I', 'M', 'M', 'R', 'E', 'A', 'R', 'Q', 'W', 'F', 'V'],
    ['K', 'B', 'K', 'D', 'J', 'F', 'C', 'C', 'X', 'Q', 'P', 'M', 'E', 'O'],
    ['H', 'D', 'Z', 'K', 'T', 'Q', 'K', 'S', 'T', 'L', 'W', 'H', 'D', 'Z']
]

result = find_largest_square(matrix)
print(result)
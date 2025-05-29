def find_largest_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    max_size = 0
    bottom_right_corner = (-1, -1)

    for i in range(rows):
        for j in range(cols):
            # Try to form a square with top-left corner at (i, j)
            for size in range(3, min(rows - i, cols - j) + 1):
                # Check if a square of this size can be formed
                char = matrix[i][j]
                valid_square = True
                # Check top and bottom rows
                for k in range(size):
                    if matrix[i][j + k] != char or matrix[i + size - 1][j + k] != char:
                        valid_square = False
                        break
                # Check left and right columns
                for k in range(size):
                    if matrix[i + k][j] != char or matrix[i + k][j + size - 1] != char:
                        valid_square = False
                        break
                if valid_square and size > max_size:
                    max_size = size
                    bottom_right_corner = (i + size - 1, j + size - 1)

    return bottom_right_corner

matrix = [
    ['S', 'M', 'Q', 'A', 'J', 'N', 'Z', 'I', 'B', 'F', 'H', 'W', 'M', 'W'],
    ['A', 'V', 'O', 'A', 'Z', 'Z', 'Z', 'M', 'X', 'B', 'C', 'N', 'M', 'V'],
    ['X', 'Z', 'Z', 'D', 'Z', 'Z', 'Z', 'S', 'V', 'J', 'N', 'R', 'D', 'L'],
    ['S', 'Z', 'T', 'O', 'Z', 'Z', 'Z', 'P', 'M', 'K', 'K', 'K', 'D', 'T'],
    ['B', 'W', 'O', 'Q', 'E', 'I', 'B', 'A', 'K', 'A', 'Q', 'X', 'N', 'T'],
    ['K', 'O', 'V', 'B', 'U', 'G', 'O', 'B', 'Q', 'S', 'F', 'R', 'U', 'J']
]

result = find_largest_square(matrix)
print(result)
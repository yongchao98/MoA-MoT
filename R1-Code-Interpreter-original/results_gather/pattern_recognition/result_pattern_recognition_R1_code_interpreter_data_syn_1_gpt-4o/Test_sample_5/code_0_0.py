def find_largest_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    max_size = 0
    bottom_right_corner = (-1, -1)

    for i in range(rows):
        for j in range(cols):
            # Try to form a square with (i, j) as the top-left corner
            for size in range(3, min(rows - i, cols - j) + 1):
                # Check if all characters on the boundary of the square are the same
                char = matrix[i][j]
                valid_square = True
                for k in range(size):
                    if matrix[i + k][j] != char or matrix[i][j + k] != char or matrix[i + size - 1][j + k] != char or matrix[i + k][j + size - 1] != char:
                        valid_square = False
                        break
                if valid_square and size > max_size:
                    max_size = size
                    bottom_right_corner = (i + size - 1, j + size - 1)

    return bottom_right_corner

matrix = [
    ['G', 'Y', 'V', 'W', 'J', 'I', 'N', 'D', 'X', 'X', 'D'],
    ['U', 'E', 'T', 'X', 'Y', 'D', 'K', 'C', 'T', 'O', 'I'],
    ['E', 'G', 'Y', 'A', 'D', 'E', 'Y', 'Y', 'Y', 'Y', 'R'],
    ['Q', 'B', 'N', 'T', 'O', 'O', 'Y', 'Y', 'Y', 'Y', 'I'],
    ['U', 'K', 'H', 'U', 'Y', 'Q', 'Y', 'Y', 'Y', 'Y', 'K'],
    ['T', 'L', 'I', 'A', 'T', 'M', 'Y', 'Y', 'Y', 'Y', 'Y'],
    ['J', 'A', 'W', 'C', 'I', 'Z', 'N', 'U', 'D', 'T', 'G'],
    ['K', 'L', 'R', 'U', 'K', 'L', 'A', 'Z', 'Y', 'C', 'K'],
    ['Z', 'S', 'J', 'R', 'T', 'R', 'X', 'V', 'M', 'Z', 'Q'],
    ['M', 'W', 'W', 'A', 'C', 'Y', 'S', 'V', 'G', 'W', 'T']
]

result = find_largest_square(matrix)
print(result)
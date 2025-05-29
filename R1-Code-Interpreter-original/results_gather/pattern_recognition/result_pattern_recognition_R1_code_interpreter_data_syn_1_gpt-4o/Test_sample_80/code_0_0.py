def find_largest_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    max_size = 0
    bottom_right_corner = (-1, -1)

    for i in range(rows):
        for j in range(cols):
            # Try to form a square with top-left corner at (i, j)
            current_char = matrix[i][j]
            for size in range(3, min(rows - i, cols - j) + 1):
                # Check if a square of this size can be formed
                valid_square = True
                # Check the right and bottom boundaries of the square
                for k in range(size):
                    if matrix[i + size - 1][j + k] != current_char or matrix[i + k][j + size - 1] != current_char:
                        valid_square = False
                        break
                if valid_square:
                    if size > max_size:
                        max_size = size
                        bottom_right_corner = (i + size - 1, j + size - 1)
                else:
                    break

    return bottom_right_corner

matrix = [
    ['Q', 'Y', 'Q', 'R', 'K', 'V', 'A', 'Q', 'Z', 'C', 'J'],
    ['S', 'K', 'F', 'C', 'N', 'U', 'X', 'O', 'X', 'X', 'I'],
    ['M', 'O', 'B', 'Q', 'Q', 'P', 'G', 'B', 'I', 'Y', 'W'],
    ['I', 'Y', 'J', 'B', 'G', 'P', 'B', 'D', 'Z', 'G', 'G'],
    ['I', 'Z', 'U', 'L', 'Y', 'P', 'M', 'M', 'M', 'M', 'T'],
    ['C', 'V', 'P', 'H', 'M', 'V', 'M', 'M', 'M', 'M', 'P'],
    ['M', 'H', 'I', 'T', 'Q', 'U', 'M', 'M', 'M', 'M', 'U'],
    ['Z', 'D', 'L', 'X', 'G', 'D', 'M', 'M', 'M', 'M', 'D'],
    ['H', 'Q', 'W', 'W', 'P', 'N', 'U', 'X', 'I', 'L', 'W'],
    ['L', 'G', 'I', 'Q', 'Y', 'P', 'E', 'W', 'B', 'U', 'Y'],
    ['B', 'Z', 'F', 'S', 'K', 'E', 'V', 'K', 'D', 'T', 'M'],
    ['B', 'Z', 'P', 'V', 'J', 'U', 'U', 'Y', 'H', 'V', 'Q']
]

result = find_largest_square(matrix)
print(result)
def find_largest_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    max_size = 0
    bottom_right_corner = (-1, -1)

    for i in range(rows):
        for j in range(cols):
            # Try to form a square with (i, j) as the top-left corner
            for size in range(3, min(rows - i, cols - j) + 1):
                # Check if a square of this size can be formed
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
                else:
                    break

    return bottom_right_corner

matrix = [
    ['U', 'D', 'L', 'G', 'V', 'L', 'V'],
    ['A', 'K', 'L', 'K', 'E', 'N', 'N'],
    ['L', 'Z', 'O', 'A', 'V', 'A', 'M'],
    ['M', 'M', 'M', 'Z', 'K', 'W', 'S'],
    ['M', 'M', 'M', 'J', 'C', 'M', 'Y'],
    ['M', 'M', 'M', 'C', 'W', 'Q', 'D'],
    ['M', 'Y', 'Q', 'P', 'I', 'X', 'Z'],
    ['A', 'F', 'F', 'U', 'D', 'U', 'G'],
    ['R', 'L', 'N', 'P', 'P', 'U', 'K'],
    ['P', 'C', 'M', 'C', 'R', 'N', 'E'],
    ['V', 'E', 'O', 'J', 'B', 'R', 'I'],
    ['X', 'Y', 'C', 'R', 'N', 'R', 'R'],
    ['U', 'C', 'F', 'Z', 'Q', 'T', 'A'],
    ['Z', 'W', 'H', 'X', 'P', 'M', 'P']
]

result = find_largest_square(matrix)
print(result)
def find_largest_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    max_size = 0
    bottom_right_corner = None

    for i in range(rows):
        for j in range(cols):
            # Try to form a square with (i, j) as the top-left corner
            for size in range(3, min(rows - i, cols - j) + 1):
                # Check if a square of this size can be formed
                char = matrix[i][j]
                valid_square = True
                # Check the right and bottom boundaries
                for k in range(size):
                    if matrix[i + k][j + size - 1] != char or matrix[i + size - 1][j + k] != char:
                        valid_square = False
                        break
                if valid_square:
                    # Update the largest square found
                    if size > max_size:
                        max_size = size
                        bottom_right_corner = (i + size - 1, j + size - 1)

    return bottom_right_corner

matrix = [
    ['M', 'N', 'V', 'J', 'R', 'J', 'G', 'U'],
    ['A', 'X', 'R', 'D', 'B', 'W', 'L', 'Z'],
    ['K', 'F', 'T', 'V', 'B', 'C', 'S', 'J'],
    ['Q', 'I', 'I', 'I', 'F', 'X', 'H', 'S'],
    ['S', 'I', 'I', 'I', 'A', 'A', 'M', 'A'],
    ['H', 'I', 'I', 'I', 'S', 'M', 'Y', 'F'],
    ['Q', 'G', 'H', 'Y', 'V', 'U', 'Y', 'H'],
    ['A', 'C', 'O', 'U', 'L', 'Z', 'J', 'C'],
    ['M', 'O', 'O', 'Q', 'N', 'N', 'O', 'H'],
    ['X', 'Y', 'C', 'H', 'H', 'U', 'E', 'R']
]

result = find_largest_square(matrix)
print(result)
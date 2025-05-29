def find_largest_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    max_size = 0
    bottom_right_corner = (-1, -1)

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
                else:
                    break

    return bottom_right_corner

matrix = [
    ['O', 'P', 'R', 'F', 'A', 'L', 'M', 'G'],
    ['C', 'Y', 'D', 'Z', 'F', 'V', 'G', 'V'],
    ['Q', 'I', 'J', 'B', 'I', 'V', 'B', 'E'],
    ['F', 'A', 'F', 'V', 'V', 'V', 'V', 'L'],
    ['D', 'B', 'Z', 'V', 'V', 'V', 'V', 'P'],
    ['Y', 'I', 'O', 'V', 'V', 'V', 'V', 'T'],
    ['K', 'K', 'J', 'V', 'V', 'V', 'V', 'Q'],
    ['C', 'V', 'O', 'F', 'E', 'Z', 'J', 'C'],
    ['J', 'F', 'D', 'V', 'C', 'F', 'S', 'A']
]

result = find_largest_square(matrix)
print(result)
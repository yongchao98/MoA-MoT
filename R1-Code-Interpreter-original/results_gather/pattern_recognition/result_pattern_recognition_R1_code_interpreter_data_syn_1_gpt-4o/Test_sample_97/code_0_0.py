def find_largest_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    max_size = 0
    bottom_right_corner = None

    # Iterate over each possible top-left corner of a square
    for i in range(rows):
        for j in range(cols):
            # Try to form a square with the top-left corner at (i, j)
            current_char = matrix[i][j]
            size = 2  # Start with size 3x3, so we need to check size 2 (0-indexed)
            while i + size < rows and j + size < cols:
                # Check if the square of current size is valid
                valid_square = True
                for k in range(size + 1):
                    if matrix[i + size][j + k] != current_char or matrix[i + k][j + size] != current_char:
                        valid_square = False
                        break
                if valid_square:
                    if size + 1 > max_size:
                        max_size = size + 1
                        bottom_right_corner = (i + size, j + size)
                else:
                    break
                size += 1

    return bottom_right_corner

matrix = [
    ['Y', 'L', 'I', 'F', 'T', 'Q', 'D', 'Y', 'F', 'I', 'V', 'X', 'L', 'V'],
    ['F', 'S', 'H', 'W', 'Y', 'J', 'X', 'R', 'I', 'T', 'I', 'R', 'B', 'H'],
    ['T', 'T', 'X', 'X', 'R', 'P', 'G', 'N', 'I', 'F', 'Q', 'T', 'N', 'A'],
    ['D', 'T', 'C', 'N', 'U', 'V', 'G', 'E', 'D', 'E', 'J', 'Z', 'M', 'R'],
    ['K', 'S', 'B', 'R', 'S', 'N', 'V', 'W', 'S', 'D', 'K', 'Y', 'W', 'J'],
    ['M', 'H', 'P', 'E', 'A', 'D', 'M', 'V', 'X', 'N', 'E', 'E', 'E', 'P'],
    ['D', 'Y', 'K', 'V', 'Q', 'H', 'N', 'U', 'A', 'X', 'E', 'E', 'E', 'A'],
    ['Z', 'G', 'D', 'R', 'D', 'M', 'F', 'Z', 'M', 'B', 'E', 'E', 'E', 'W']
]

result = find_largest_square(matrix)
print(result)
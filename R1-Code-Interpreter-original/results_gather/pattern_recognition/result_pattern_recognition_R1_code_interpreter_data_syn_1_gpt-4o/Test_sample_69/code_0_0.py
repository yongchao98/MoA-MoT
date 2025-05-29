def find_largest_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    max_size = 0
    bottom_right_corner = (-1, -1)

    for i in range(rows):
        for j in range(cols):
            # Try to form a square with top-left corner at (i, j)
            current_char = matrix[i][j]
            size = 2  # Start checking from size 3x3
            while i + size < rows and j + size < cols:
                valid_square = True
                # Check the new row and column added to the square
                for k in range(size + 1):
                    if matrix[i + size][j + k] != current_char or matrix[i + k][j + size] != current_char:
                        valid_square = False
                        break
                if valid_square:
                    size += 1
                else:
                    break

            # Update the largest square found
            if size > max_size:
                max_size = size
                bottom_right_corner = (i + size - 1, j + size - 1)

    return bottom_right_corner

matrix = [
    ['X', 'A', 'B', 'C', 'F', 'F', 'F'],
    ['E', 'T', 'Y', 'G', 'F', 'F', 'F'],
    ['U', 'G', 'P', 'S', 'F', 'F', 'F'],
    ['W', 'E', 'A', 'V', 'P', 'A', 'I'],
    ['J', 'W', 'E', 'F', 'A', 'Q', 'V'],
    ['F', 'S', 'G', 'J', 'K', 'Y', 'Y'],
    ['S', 'P', 'Q', 'N', 'B', 'U', 'K'],
    ['W', 'A', 'G', 'Z', 'Z', 'H', 'Z']
]

result = find_largest_square(matrix)
print(result)
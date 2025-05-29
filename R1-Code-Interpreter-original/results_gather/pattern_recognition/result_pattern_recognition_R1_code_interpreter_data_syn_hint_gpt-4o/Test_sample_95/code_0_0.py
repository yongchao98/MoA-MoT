def find_largest_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    max_size = 0
    bottom_right_corner = None

    for i in range(rows):
        for j in range(cols):
            # Try to form a square with top-left corner at (i, j)
            current_char = matrix[i][j]
            size = 2  # Start checking from 3x3 squares
            while i + size < rows and j + size < cols:
                # Check if the square of size (size+1)x(size+1) is valid
                valid_square = True
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
    ['U', 'Y', 'R', 'Z', 'D', 'D', 'D', 'E'],
    ['H', 'S', 'Z', 'Z', 'D', 'D', 'D', 'C'],
    ['H', 'H', 'O', 'Z', 'D', 'D', 'D', 'U'],
    ['C', 'O', 'C', 'B', 'E', 'O', 'Q', 'D'],
    ['S', 'W', 'F', 'W', 'N', 'R', 'R', 'G'],
    ['Y', 'V', 'A', 'W', 'O', 'J', 'D', 'W'],
    ['L', 'N', 'U', 'I', 'T', 'G', 'A', 'Z']
]

result = find_largest_square(matrix)
print(result)
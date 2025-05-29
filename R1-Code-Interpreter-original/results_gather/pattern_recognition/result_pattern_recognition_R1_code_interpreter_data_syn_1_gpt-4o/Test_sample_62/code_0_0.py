def find_largest_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    max_size = 0
    bottom_right_corner = (-1, -1)

    # Iterate over each possible top-left corner of a square
    for i in range(rows):
        for j in range(cols):
            # Try to form squares of increasing size
            for size in range(3, min(rows - i, cols - j) + 1):
                # Check if all characters in the square are the same
                char = matrix[i][j]
                is_square = True
                for x in range(i, i + size):
                    for y in range(j, j + size):
                        if matrix[x][y] != char:
                            is_square = False
                            break
                    if not is_square:
                        break
                
                # If a valid square is found, update the largest square
                if is_square and size > max_size:
                    max_size = size
                    bottom_right_corner = (i + size - 1, j + size - 1)

    return bottom_right_corner

# Define the matrix
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

# Find the bottom-right corner of the largest square
result = find_largest_square(matrix)
print(result)
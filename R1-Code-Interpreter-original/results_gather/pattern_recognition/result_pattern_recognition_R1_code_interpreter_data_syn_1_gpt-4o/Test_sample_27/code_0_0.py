def find_largest_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    max_square_size = 0
    bottom_right_corner = None

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
                
                # If a valid square is found, update the largest square info
                if is_square and size > max_square_size:
                    max_square_size = size
                    bottom_right_corner = (i + size - 1, j + size - 1)

    return bottom_right_corner

# Define the matrix
matrix = [
    ['D', 'N', 'E', 'L', 'R', 'Z', 'R', 'R', 'R', 'O', 'P', 'Z', 'R'],
    ['A', 'T', 'Q', 'D', 'M', 'M', 'R', 'R', 'R', 'H', 'D', 'P', 'K'],
    ['U', 'U', 'U', 'A', 'V', 'D', 'R', 'R', 'R', 'Q', 'T', 'F', 'B'],
    ['T', 'M', 'D', 'F', 'O', 'W', 'Y', 'C', 'O', 'N', 'P', 'G', 'D'],
    ['L', 'V', 'S', 'N', 'L', 'X', 'X', 'B', 'Y', 'J', 'A', 'L', 'F'],
    ['X', 'D', 'O', 'B', 'W', 'J', 'G', 'H', 'V', 'U', 'A', 'M', 'S'],
    ['G', 'N', 'M', 'U', 'G', 'Z', 'P', 'R', 'U', 'D', 'J', 'B', 'Z'],
    ['Z', 'G', 'M', 'A', 'Y', 'F', 'D', 'J', 'Z', 'L', 'O', 'F', 'N'],
    ['C', 'Q', 'L', 'E', 'U', 'A', 'T', 'N', 'Q', 'S', 'V', 'O', 'O'],
    ['L', 'A', 'X', 'M', 'X', 'D', 'R', 'L', 'X', 'C', 'T', 'A', 'S'],
    ['V', 'W', 'F', 'K', 'W', 'E', 'P', 'K', 'V', 'W', 'D', 'U', 'Z'],
    ['C', 'R', 'A', 'F', 'R', 'R', 'S', 'N', 'O', 'Q', 'C', 'G', 'S']
]

# Find the bottom-right corner of the largest square
result = find_largest_square(matrix)
print(result)
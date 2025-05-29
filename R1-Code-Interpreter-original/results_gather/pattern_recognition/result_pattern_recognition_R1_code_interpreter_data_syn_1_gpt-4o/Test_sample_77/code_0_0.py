def find_largest_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    max_size = 0
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
                if is_square and size > max_size:
                    max_size = size
                    bottom_right_corner = (i + size - 1, j + size - 1)

    return bottom_right_corner

# Define the matrix
matrix = [
    ['M', 'U', 'Z', 'H', 'H', 'H', 'H', 'N', 'T'],
    ['A', 'H', 'O', 'H', 'H', 'H', 'H', 'M', 'N'],
    ['R', 'A', 'J', 'H', 'H', 'H', 'H', 'Q', 'G'],
    ['K', 'J', 'X', 'H', 'H', 'H', 'H', 'V', 'L'],
    ['V', 'Y', 'L', 'Q', 'E', 'K', 'D', 'P', 'Y'],
    ['O', 'J', 'P', 'E', 'X', 'I', 'P', 'V', 'R'],
    ['T', 'H', 'V', 'B', 'Z', 'V', 'P', 'A', 'A'],
    ['X', 'D', 'C', 'H', 'S', 'D', 'H', 'S', 'C'],
    ['P', 'R', 'Q', 'Y', 'Q', 'O', 'Y', 'N', 'C'],
    ['X', 'A', 'O', 'M', 'R', 'K', 'D', 'X', 'I'],
    ['C', 'A', 'H', 'X', 'G', 'L', 'S', 'L', 'F'],
    ['Q', 'O', 'R', 'U', 'P', 'N', 'B', 'I', 'X'],
    ['D', 'R', 'Y', 'S', 'W', 'Q', 'W', 'W', 'E'],
    ['Q', 'V', 'Y', 'Y', 'A', 'Q', 'V', 'H', 'O'],
    ['N', 'I', 'V', 'N', 'V', 'X', 'O', 'Q', 'D']
]

# Find the bottom-right corner of the largest square
result = find_largest_square(matrix)
print(result)
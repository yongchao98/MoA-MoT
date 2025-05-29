def find_largest_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    max_size = 0
    bottom_right_corner = (-1, -1)

    # Iterate over each cell as a potential top-left corner of a square
    for i in range(rows):
        for j in range(cols):
            # Check for squares of increasing size
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
    ['L', 'D', 'H', 'X', 'Q', 'G', 'D', 'V', 'E', 'S'],
    ['E', 'H', 'Y', 'T', 'D', 'Y', 'R', 'T', 'N', 'M'],
    ['K', 'K', 'K', 'K', 'A', 'M', 'E', 'M', 'M', 'K'],
    ['K', 'K', 'K', 'K', 'E', 'W', 'X', 'D', 'L', 'N'],
    ['K', 'K', 'K', 'K', 'Z', 'D', 'U', 'D', 'T', 'N'],
    ['K', 'K', 'K', 'K', 'Y', 'N', 'V', 'J', 'V', 'F'],
    ['Z', 'C', 'T', 'J', 'M', 'G', 'V', 'D', 'X', 'P'],
    ['G', 'D', 'A', 'J', 'R', 'U', 'Y', 'H', 'K', 'L'],
    ['S', 'U', 'N', 'E', 'R', 'A', 'Y', 'X', 'N', 'Y'],
    ['K', 'Q', 'I', 'O', 'P', 'N', 'R', 'G', 'N', 'V']
]

# Find the bottom-right corner of the largest square
result = find_largest_square(matrix)
print(result)
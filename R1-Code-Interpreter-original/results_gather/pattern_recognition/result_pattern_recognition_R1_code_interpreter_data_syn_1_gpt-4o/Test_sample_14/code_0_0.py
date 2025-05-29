def find_largest_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    max_square_size = 0
    bottom_right_corner = None

    for i in range(rows):
        for j in range(cols):
            # Try to form squares with top-left corner at (i, j)
            for size in range(3, min(rows - i, cols - j) + 1):
                # Check if all characters in the square are the same
                char = matrix[i][j]
                valid_square = True
                for x in range(i, i + size):
                    for y in range(j, j + size):
                        if matrix[x][y] != char:
                            valid_square = False
                            break
                    if not valid_square:
                        break
                
                if valid_square:
                    if size > max_square_size:
                        max_square_size = size
                        bottom_right_corner = (i + size - 1, j + size - 1)

    return bottom_right_corner

matrix = [
    ['G', 'H', 'Y', 'J', 'I', 'I', 'J', 'D', 'S', 'L', 'E'],
    ['T', 'F', 'G', 'S', 'U', 'G', 'I', 'D', 'W', 'V', 'K'],
    ['P', 'O', 'D', 'A', 'Y', 'P', 'F', 'X', 'C', 'Q', 'X'],
    ['Y', 'U', 'H', 'M', 'Y', 'D', 'D', 'D', 'D', 'D', 'M'],
    ['B', 'D', 'H', 'D', 'P', 'D', 'D', 'D', 'D', 'D', 'A'],
    ['E', 'I', 'N', 'O', 'N', 'D', 'D', 'D', 'D', 'D', 'I'],
    ['U', 'J', 'C', 'U', 'B', 'D', 'D', 'D', 'D', 'D', 'V'],
    ['Y', 'P', 'K', 'I', 'D', 'D', 'D', 'D', 'D', 'D', 'K'],
    ['E', 'U', 'G', 'U', 'D', 'J', 'F', 'O', 'N', 'G', 'L'],
    ['P', 'F', 'F', 'O', 'V', 'R', 'F', 'M', 'Y', 'O', 'D'],
    ['G', 'M', 'X', 'Z', 'E', 'D', 'J', 'B', 'U', 'W', 'J'],
    ['B', 'V', 'A', 'A', 'M', 'C', 'K', 'P', 'Z', 'S', 'U'],
    ['G', 'A', 'A', 'G', 'X', 'J', 'T', 'L', 'R', 'B', 'W'],
    ['I', 'N', 'Q', 'I', 'U', 'L', 'S', 'L', 'F', 'E', 'D']
]

result = find_largest_square(matrix)
print(result)
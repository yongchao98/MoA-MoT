def find_largest_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    max_size = 0
    bottom_right_corner = (-1, -1)

    for i in range(rows):
        for j in range(cols):
            # Try to form a square with (i, j) as the top-left corner
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
                
                if is_square and size > max_size:
                    max_size = size
                    bottom_right_corner = (i + size - 1, j + size - 1)

    return bottom_right_corner

matrix = [
    ['C', 'S', 'E', 'G', 'Y', 'Q', 'J', 'Y', 'W', 'R', 'Y'],
    ['O', 'F', 'D', 'Z', 'Q', 'H', 'D', 'I', 'W', 'Q', 'O'],
    ['T', 'F', 'F', 'F', 'F', 'S', 'U', 'U', 'Y', 'A', 'K'],
    ['U', 'F', 'F', 'F', 'F', 'B', 'Y', 'X', 'W', 'C', 'V'],
    ['C', 'F', 'F', 'F', 'F', 'P', 'S', 'P', 'O', 'O', 'Z'],
    ['S', 'F', 'F', 'F', 'F', 'X', 'L', 'F', 'T', 'J', 'B'],
    ['G', 'A', 'P', 'E', 'U', 'K', 'L', 'L', 'K', 'K', 'W'],
    ['X', 'W', 'G', 'O', 'N', 'S', 'Y', 'Y', 'E', 'E', 'Y'],
    ['E', 'L', 'R', 'Z', 'P', 'J', 'U', 'U', 'V', 'P', 'C'],
    ['H', 'K', 'X', 'Z', 'Y', 'O', 'N', 'U', 'K', 'F', 'F'],
    ['A', 'X', 'J', 'K', 'K', 'D', 'O', 'R', 'W', 'P', 'Y'],
    ['P', 'M', 'R', 'F', 'X', 'U', 'J', 'Z', 'W', 'F', 'P'],
    ['R', 'J', 'R', 'J', 'P', 'L', 'A', 'U', 'Z', 'S', 'Y'],
    ['T', 'V', 'E', 'F', 'Z', 'J', 'H', 'M', 'Z', 'W', 'Y']
]

result = find_largest_square(matrix)
print(result)
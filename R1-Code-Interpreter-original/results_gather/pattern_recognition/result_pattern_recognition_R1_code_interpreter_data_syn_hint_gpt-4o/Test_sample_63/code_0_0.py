def find_largest_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    max_square_size = 0
    bottom_right_corner = None

    for i in range(rows):
        for j in range(cols):
            # Check for squares with top-left corner at (i, j)
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
                if is_square and size > max_square_size:
                    max_square_size = size
                    bottom_right_corner = (i + size - 1, j + size - 1)

    return bottom_right_corner

matrix = [
    ['O', 'J', 'N', 'F', 'Z', 'X', 'B', 'U', 'R'],
    ['I', 'W', 'Y', 'B', 'N', 'U', 'O', 'S', 'Q'],
    ['O', 'I', 'B', 'K', 'R', 'U', 'Q', 'H', 'S'],
    ['Z', 'Z', 'H', 'U', 'Y', 'Q', 'S', 'R', 'A'],
    ['C', 'C', 'K', 'N', 'T', 'B', 'B', 'C', 'V'],
    ['T', 'Y', 'W', 'S', 'N', 'U', 'I', 'W', 'E'],
    ['Q', 'O', 'O', 'V', 'X', 'B', 'N', 'V', 'T'],
    ['X', 'H', 'F', 'N', 'N', 'N', 'P', 'P', 'J'],
    ['V', 'C', 'H', 'N', 'N', 'N', 'A', 'O', 'C'],
    ['H', 'B', 'E', 'N', 'N', 'N', 'G', 'Y', 'H'],
    ['Q', 'K', 'I', 'E', 'T', 'J', 'D', 'E', 'F'],
    ['U', 'C', 'S', 'A', 'A', 'Q', 'P', 'L', 'U'],
    ['L', 'O', 'E', 'R', 'O', 'E', 'H', 'C', 'F'],
    ['E', 'R', 'S', 'V', 'Q', 'I', 'D', 'C', 'G']
]

result = find_largest_square(matrix)
print(result)
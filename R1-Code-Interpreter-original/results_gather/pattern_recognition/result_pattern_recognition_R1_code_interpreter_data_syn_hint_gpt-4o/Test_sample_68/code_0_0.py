def find_largest_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    max_size = 0
    bottom_right_corner = None

    for i in range(rows):
        for j in range(cols):
            # Try to form a square with top-left corner at (i, j)
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
    ['T', 'X', 'A', 'L', 'J', 'V', 'X', 'Y', 'Y', 'L', 'R', 'U', 'X', 'I'],
    ['E', 'G', 'D', 'U', 'Y', 'T', 'C', 'F', 'F', 'F', 'F', 'L', 'I', 'U'],
    ['K', 'X', 'S', 'B', 'E', 'R', 'D', 'F', 'F', 'F', 'F', 'U', 'S', 'R'],
    ['F', 'G', 'U', 'K', 'P', 'F', 'O', 'F', 'F', 'F', 'F', 'F', 'C', 'E'],
    ['C', 'K', 'W', 'Z', 'C', 'V', 'Y', 'F', 'F', 'F', 'F', 'F', 'K', 'D'],
    ['H', 'E', 'N', 'J', 'Y', 'T', 'O', 'V', 'M', 'U', 'O', 'Y', 'J', 'J'],
    ['R', 'A', 'Y', 'Y', 'N', 'K', 'K', 'I', 'K', 'I', 'O', 'P', 'D', 'P'],
    ['E', 'M', 'G', 'B', 'E', 'E', 'G', 'X', 'G', 'P', 'T', 'R', 'K', 'K']
]

result = find_largest_square(matrix)
print(result)
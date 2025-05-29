def find_largest_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    max_size = 0
    bottom_right_corner = None

    for i in range(rows):
        for j in range(cols):
            # Try to form a square with (i, j) as the top-left corner
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
                
                if valid_square and size > max_size:
                    max_size = size
                    bottom_right_corner = (i + size - 1, j + size - 1)

    return bottom_right_corner

matrix = [
    ['G', 'Y', 'A', 'D', 'H', 'G', 'C', 'P', 'H', 'F', 'E', 'E', 'F', 'G', 'F'],
    ['X', 'K', 'G', 'F', 'T', 'E', 'D', 'D', 'E', 'Y', 'H', 'N', 'J', 'R', 'P'],
    ['V', 'L', 'B', 'F', 'A', 'I', 'I', 'I', 'I', 'H', 'O', 'N', 'O', 'R', 'H'],
    ['E', 'N', 'G', 'U', 'Z', 'I', 'I', 'I', 'I', 'W', 'J', 'O', 'M', 'V', 'F'],
    ['J', 'E', 'M', 'F', 'U', 'I', 'I', 'I', 'I', 'A', 'Y', 'B', 'T', 'F', 'E'],
    ['D', 'R', 'A', 'H', 'R', 'I', 'I', 'I', 'I', 'X', 'R', 'M', 'I', 'S', 'R'],
    ['O', 'O', 'O', 'C', 'X', 'X', 'O', 'R', 'F', 'F', 'K', 'L', 'Y', 'G', 'N'],
    ['P', 'R', 'X', 'Y', 'M', 'G', 'W', 'L', 'G', 'A', 'R', 'B', 'G', 'I', 'F']
]

result = find_largest_square(matrix)
print(result)
def find_largest_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    max_size = 0
    bottom_right_corner = (-1, -1)

    for i in range(rows):
        for j in range(cols):
            # Try to form a square with top-left corner at (i, j)
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
    ['M', 'Y', 'E', 'O', 'Q', 'B', 'O', 'N'],
    ['P', 'P', 'D', 'M', 'B', 'F', 'X', 'P'],
    ['J', 'E', 'X', 'O', 'O', 'O', 'O', 'O'],
    ['V', 'G', 'C', 'O', 'O', 'O', 'O', 'W'],
    ['H', 'Z', 'K', 'O', 'O', 'O', 'O', 'O'],
    ['J', 'V', 'L', 'O', 'O', 'O', 'O', 'E'],
    ['K', 'Y', 'N', 'Z', 'Q', 'M', 'F', 'B'],
    ['F', 'R', 'T', 'N', 'W', 'X', 'O', 'T'],
    ['T', 'A', 'D', 'X', 'U', 'N', 'U', 'U'],
    ['X', 'B', 'X', 'E', 'Y', 'M', 'R', 'W']
]

result = find_largest_square(matrix)
print(result)
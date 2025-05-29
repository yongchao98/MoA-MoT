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
    ['O', 'C', 'D', 'G', 'G', 'Q', 'J', 'Q', 'U', 'K'],
    ['X', 'U', 'E', 'C', 'B', 'R', 'V', 'S', 'K', 'M'],
    ['E', 'G', 'Z', 'R', 'R', 'J', 'R', 'Y', 'O', 'T'],
    ['W', 'U', 'J', 'J', 'V', 'P', 'G', 'C', 'Z', 'X'],
    ['K', 'T', 'F', 'L', 'P', 'P', 'N', 'H', 'N', 'I'],
    ['K', 'K', 'K', 'Z', 'E', 'Q', 'E', 'Q', 'C', 'Y'],
    ['K', 'K', 'K', 'U', 'U', 'W', 'U', 'K', 'S', 'B'],
    ['K', 'K', 'K', 'O', 'B', 'O', 'B', 'F', 'S', 'V'],
    ['B', 'O', 'L', 'B', 'L', 'Z', 'T', 'Z', 'F', 'T'],
    ['N', 'B', 'T', 'X', 'O', 'Z', 'N', 'L', 'L', 'Y']
]

result = find_largest_square(matrix)
print(result)
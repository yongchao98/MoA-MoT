def find_largest_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    max_size = 0
    bottom_right_corner = None

    for i in range(rows):
        for j in range(cols):
            # Check for squares starting from size 3
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
    ['U', 'O', 'W', 'N', 'E', 'P', 'V', 'H'],
    ['T', 'T', 'T', 'T', 'T', 'N', 'N', 'C'],
    ['T', 'T', 'T', 'B', 'K', 'P', 'I', 'Q'],
    ['T', 'T', 'T', 'Y', 'T', 'K', 'H', 'C'],
    ['O', 'U', 'Y', 'X', 'E', 'K', 'G', 'E'],
    ['E', 'K', 'O', 'E', 'X', 'F', 'B', 'X'],
    ['P', 'X', 'T', 'K', 'L', 'W', 'H', 'U']
]

result = find_largest_square(matrix)
print(result)
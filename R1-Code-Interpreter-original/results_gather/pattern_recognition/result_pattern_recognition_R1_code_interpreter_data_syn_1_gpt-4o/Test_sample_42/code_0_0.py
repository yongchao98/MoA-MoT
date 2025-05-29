def find_largest_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    max_size = 0
    bottom_right_corner = (-1, -1)

    for i in range(rows):
        for j in range(cols):
            # Check for squares starting from size 3x3
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
    ['W', 'E', 'B', 'P', 'M', 'L'],
    ['E', 'R', 'Q', 'T', 'S', 'L'],
    ['Z', 'V', 'V', 'V', 'D', 'H'],
    ['R', 'V', 'V', 'V', 'F', 'A'],
    ['K', 'V', 'V', 'V', 'N', 'F'],
    ['R', 'K', 'R', 'V', 'D', 'J'],
    ['J', 'W', 'T', 'E', 'W', 'J'],
    ['E', 'A', 'M', 'R', 'E', 'K']
]

result = find_largest_square(matrix)
print(result)
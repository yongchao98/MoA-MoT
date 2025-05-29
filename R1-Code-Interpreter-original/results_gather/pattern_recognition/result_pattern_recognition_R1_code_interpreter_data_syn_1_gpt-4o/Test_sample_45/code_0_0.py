def find_largest_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    max_square_size = min(rows, cols)
    result = None

    for size in range(3, max_square_size + 1):
        for i in range(rows - size + 1):
            for j in range(cols - size + 1):
                char = matrix[i][j]
                is_square = True
                for x in range(size):
                    for y in range(size):
                        if matrix[i + x][j + y] != char:
                            is_square = False
                            break
                    if not is_square:
                        break
                if is_square:
                    result = (i + size - 1, j + size - 1)

    return result

matrix = [
    ['Y', 'Y', 'J', 'F', 'J', 'V', 'M'],
    ['C', 'O', 'D', 'K', 'W', 'M', 'P'],
    ['F', 'B', 'K', 'K', 'K', 'A', 'O'],
    ['U', 'P', 'K', 'K', 'K', 'X', 'H'],
    ['F', 'R', 'K', 'K', 'K', 'J', 'Y'],
    ['S', 'Y', 'M', 'C', 'V', 'P', 'G']
]

result = find_largest_square(matrix)
print(result)
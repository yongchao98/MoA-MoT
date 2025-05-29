def find_largest_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    max_size = min(rows, cols)
    result = None

    for size in range(3, max_size + 1):
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
    ['S', 'L', 'R', 'S', 'D', 'A', 'J', 'D', 'K', 'K'],
    ['D', 'Z', 'O', 'T', 'M', 'N', 'O', 'I', 'I', 'L'],
    ['T', 'J', 'X', 'K', 'F', 'F', 'F', 'F', 'F', 'W'],
    ['Q', 'Z', 'A', 'A', 'F', 'F', 'F', 'F', 'F', 'Q'],
    ['L', 'W', 'R', 'U', 'F', 'F', 'F', 'F', 'F', 'S'],
    ['Y', 'M', 'X', 'M', 'F', 'F', 'F', 'F', 'F', 'N'],
    ['P', 'E', 'I', 'E', 'F', 'F', 'F', 'F', 'F', 'F'],
    ['X', 'F', 'K', 'C', 'W', 'J', 'R', 'W', 'N', 'P'],
    ['S', 'B', 'P', 'M', 'V', 'J', 'S', 'S', 'L', 'L'],
    ['C', 'S', 'M', 'U', 'Z', 'U', 'U', 'T', 'F', 'W'],
    ['E', 'R', 'J', 'A', 'S', 'P', 'T', 'Z', 'L', 'S'],
    ['X', 'C', 'L', 'R', 'Y', 'L', 'B', 'X', 'W', 'I']
]

result = find_largest_square(matrix)
print(result)
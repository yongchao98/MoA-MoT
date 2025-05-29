def find_largest_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    max_size = min(rows, cols)
    largest_square = None

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
                    largest_square = (i + size - 1, j + size - 1)

    return largest_square

matrix = [
    ['E', 'N', 'W', 'Y', 'F', 'A', 'R'],
    ['O', 'C', 'Z', 'U', 'J', 'D', 'Y'],
    ['Z', 'C', 'C', 'C', 'C', 'C', 'P'],
    ['N', 'J', 'S', 'C', 'C', 'C', 'A'],
    ['M', 'H', 'H', 'C', 'C', 'C', 'S'],
    ['P', 'W', 'A', 'A', 'P', 'M', 'U'],
    ['R', 'P', 'S', 'F', 'M', 'P', 'B'],
    ['H', 'T', 'B', 'C', 'X', 'J', 'B'],
    ['O', 'Z', 'X', 'A', 'P', 'E', 'V']
]

result = find_largest_square(matrix)
print(result)
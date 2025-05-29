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
    ['G', 'H', 'Y', 'J', 'I', 'I', 'J', 'D', 'S', 'L', 'E'],
    ['T', 'F', 'G', 'S', 'U', 'G', 'I', 'D', 'W', 'V', 'K'],
    ['P', 'O', 'D', 'A', 'Y', 'P', 'F', 'X', 'C', 'Q', 'X'],
    ['Y', 'U', 'H', 'M', 'Y', 'D', 'D', 'D', 'D', 'D', 'M'],
    ['B', 'D', 'H', 'D', 'P', 'D', 'D', 'D', 'D', 'D', 'A'],
    ['E', 'I', 'N', 'O', 'N', 'D', 'D', 'D', 'D', 'D', 'I'],
    ['U', 'J', 'C', 'U', 'B', 'D', 'D', 'D', 'D', 'D', 'V'],
    ['Y', 'P', 'K', 'I', 'D', 'D', 'D', 'D', 'D', 'D', 'K'],
    ['E', 'U', 'G', 'U', 'D', 'J', 'F', 'O', 'N', 'G', 'L'],
    ['P', 'F', 'F', 'O', 'V', 'R', 'F', 'M', 'Y', 'O', 'D'],
    ['G', 'M', 'X', 'Z', 'E', 'D', 'J', 'B', 'U', 'W', 'J'],
    ['B', 'V', 'A', 'A', 'M', 'C', 'K', 'P', 'Z', 'S', 'U'],
    ['G', 'A', 'A', 'G', 'X', 'J', 'T', 'L', 'R', 'B', 'W'],
    ['I', 'N', 'Q', 'I', 'U', 'L', 'S', 'L', 'F', 'E', 'D']
]

result = find_largest_square(matrix)
print(result)
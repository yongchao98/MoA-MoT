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
    ['G', 'Y', 'V', 'W', 'J', 'I', 'N', 'D', 'X', 'X', 'D'],
    ['U', 'E', 'T', 'X', 'Y', 'D', 'K', 'C', 'T', 'O', 'I'],
    ['E', 'G', 'Y', 'A', 'D', 'E', 'Y', 'Y', 'Y', 'Y', 'R'],
    ['Q', 'B', 'N', 'T', 'O', 'O', 'Y', 'Y', 'Y', 'Y', 'I'],
    ['U', 'K', 'H', 'U', 'Y', 'Q', 'Y', 'Y', 'Y', 'Y', 'K'],
    ['T', 'L', 'I', 'A', 'T', 'M', 'Y', 'Y', 'Y', 'Y', 'Y'],
    ['J', 'A', 'W', 'C', 'I', 'Z', 'N', 'U', 'D', 'T', 'G'],
    ['K', 'L', 'R', 'U', 'K', 'L', 'A', 'Z', 'Y', 'C', 'K'],
    ['Z', 'S', 'J', 'R', 'T', 'R', 'X', 'V', 'M', 'Z', 'Q'],
    ['M', 'W', 'W', 'A', 'C', 'Y', 'S', 'V', 'G', 'W', 'T']
]

result = find_largest_square(matrix)
print(result)
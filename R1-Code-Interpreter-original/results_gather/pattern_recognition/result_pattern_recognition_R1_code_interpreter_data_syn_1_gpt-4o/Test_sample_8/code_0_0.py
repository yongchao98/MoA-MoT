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
    ['G', 'I', 'V', 'Q', 'I', 'G'],
    ['U', 'V', 'L', 'F', 'V', 'W'],
    ['C', 'Y', 'V', 'V', 'J', 'T'],
    ['R', 'M', 'J', 'J', 'Q', 'R'],
    ['C', 'O', 'T', 'K', 'T', 'X'],
    ['Z', 'M', 'B', 'B', 'B', 'T'],
    ['L', 'F', 'B', 'B', 'B', 'X'],
    ['H', 'R', 'B', 'B', 'B', 'Q'],
    ['T', 'H', 'W', 'A', 'T', 'Q'],
    ['G', 'I', 'U', 'M', 'U', 'V']
]

result = find_largest_square(matrix)
print(result)
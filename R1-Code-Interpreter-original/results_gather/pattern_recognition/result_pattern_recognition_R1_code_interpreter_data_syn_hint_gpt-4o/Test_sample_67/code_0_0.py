def find_largest_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    max_size = 0
    bottom_right_corner = None

    for i in range(rows):
        for j in range(cols):
            # Check for possible square sizes starting from 3x3
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
    ['P', 'K', 'S', 'C', 'F', 'S', 'E', 'U', 'W', 'G', 'B', 'T', 'T'],
    ['L', 'B', 'N', 'L', 'L', 'O', 'A', 'G', 'Q', 'P', 'M', 'H', 'R'],
    ['F', 'O', 'S', 'D', 'P', 'H', 'I', 'X', 'V', 'N', 'X', 'M', 'F'],
    ['R', 'E', 'H', 'T', 'U', 'F', 'I', 'J', 'E', 'M', 'R', 'P', 'O'],
    ['T', 'A', 'F', 'F', 'M', 'X', 'J', 'I', 'C', 'N', 'Z', 'Y', 'H'],
    ['Q', 'D', 'U', 'M', 'X', 'D', 'M', 'N', 'T', 'N', 'H', 'V', 'Q'],
    ['H', 'U', 'M', 'D', 'Q', 'L', 'R', 'R', 'V', 'V', 'I', 'W', 'Y'],
    ['A', 'I', 'B', 'Y', 'S', 'L', 'D', 'E', 'X', 'W', 'E', 'N', 'K'],
    ['J', 'A', 'H', 'T', 'H', 'G', 'K', 'O', 'L', 'H', 'C', 'Y', 'W'],
    ['T', 'X', 'U', 'I', 'G', 'W', 'X', 'K', 'H', 'H', 'B', 'F', 'D'],
    ['L', 'A', 'U', 'I', 'V', 'V', 'P', 'T', 'Y', 'X', 'F', 'F', 'J'],
    ['F', 'S', 'R', 'R', 'R', 'K', 'G', 'Y', 'L', 'Z', 'M', 'V', 'P'],
    ['R', 'I', 'R', 'R', 'R', 'P', 'S', 'X', 'H', 'Y', 'V', 'T', 'P'],
    ['V', 'S', 'R', 'R', 'R', 'U', 'U', 'T', 'K', 'J', 'N', 'N', 'S'],
    ['I', 'X', 'M', 'P', 'K', 'Q', 'C', 'I', 'F', 'U', 'Q', 'Y', 'J']
]

result = find_largest_square(matrix)
print(result)
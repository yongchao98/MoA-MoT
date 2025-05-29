def find_largest_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    max_size = 0
    bottom_right_corner = (-1, -1)

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
    ['E', 'O', 'T', 'I', 'A', 'S'],
    ['D', 'Y', 'R', 'B', 'U', 'O'],
    ['Y', 'M', 'S', 'Z', 'Z', 'Z'],
    ['N', 'G', 'H', 'Z', 'Z', 'Z'],
    ['Q', 'S', 'X', 'Z', 'Z', 'Z'],
    ['V', 'P', 'M', 'K', 'Z', 'J'],
    ['E', 'X', 'X', 'M', 'F', 'N']
]

result = find_largest_square(matrix)
print(result)
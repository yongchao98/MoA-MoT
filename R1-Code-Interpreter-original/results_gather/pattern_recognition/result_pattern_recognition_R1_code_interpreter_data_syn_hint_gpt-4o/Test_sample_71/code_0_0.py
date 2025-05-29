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
    ['U', 'L', 'M', 'H', 'P', 'C', 'C', 'S', 'T', 'Y', 'L', 'P'],
    ['A', 'C', 'L', 'O', 'U', 'R', 'G', 'G', 'X', 'I', 'O', 'T'],
    ['H', 'F', 'W', 'S', 'C', 'B', 'C', 'W', 'D', 'V', 'U', 'G'],
    ['K', 'M', 'M', 'M', 'M', 'Q', 'S', 'Y', 'Z', 'U', 'Q', 'E'],
    ['C', 'M', 'M', 'M', 'M', 'O', 'O', 'E', 'S', 'A', 'P', 'I'],
    ['N', 'M', 'M', 'M', 'M', 'F', 'W', 'V', 'Z', 'M', 'D', 'Q'],
    ['Q', 'M', 'M', 'M', 'M', 'S', 'W', 'R', 'Y', 'T', 'O', 'N'],
    ['N', 'Z', 'C', 'Z', 'Z', 'V', 'N', 'C', 'R', 'L', 'F', 'P']
]

result = find_largest_square(matrix)
print(result)
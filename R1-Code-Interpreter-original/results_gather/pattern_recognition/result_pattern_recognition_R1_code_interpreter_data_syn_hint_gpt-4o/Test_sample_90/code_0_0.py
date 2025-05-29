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
    ['I', 'L', 'M', 'W', 'Y', 'Q', 'C', 'P', 'Z', 'X', 'J', 'M'],
    ['A', 'I', 'L', 'Q', 'Q', 'Q', 'Q', 'Q', 'X', 'H', 'V', 'B'],
    ['Q', 'I', 'G', 'Q', 'Q', 'Q', 'Q', 'Q', 'M', 'R', 'G', 'T'],
    ['D', 'K', 'V', 'Q', 'Q', 'Q', 'Q', 'Q', 'Z', 'U', 'C', 'R'],
    ['H', 'W', 'S', 'Q', 'Q', 'Q', 'Q', 'Q', 'N', 'U', 'T', 'E'],
    ['N', 'T', 'C', 'Q', 'Q', 'Q', 'Q', 'Q', 'V', 'L', 'E', 'Z'],
    ['S', 'V', 'B', 'A', 'C', 'Y', 'M', 'Q', 'S', 'N', 'P', 'E'],
    ['M', 'D', 'U', 'L', 'P', 'V', 'G', 'N', 'K', 'X', 'Y', 'P'],
    ['N', 'T', 'Y', 'D', 'R', 'D', 'V', 'K', 'U', 'L', 'K', 'F'],
    ['V', 'X', 'P', 'A', 'F', 'H', 'O', 'B', 'W', 'U', 'K', 'K'],
    ['E', 'P', 'E', 'J', 'Z', 'P', 'F', 'O', 'L', 'N', 'F', 'T'],
    ['G', 'A', 'A', 'X', 'M', 'R', 'K', 'P', 'P', 'A', 'C', 'P']
]

result = find_largest_square(matrix)
print(result)
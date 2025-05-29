def find_largest_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    max_size = 0
    bottom_right_corner = None

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
    ['R', 'L', 'L', 'L', 'Q', 'A', 'K', 'G', 'G', 'V'],
    ['H', 'L', 'L', 'L', 'R', 'U', 'Z', 'X', 'P', 'H'],
    ['N', 'L', 'L', 'L', 'N', 'J', 'F', 'V', 'P', 'J'],
    ['E', 'V', 'P', 'B', 'O', 'P', 'H', 'H', 'T', 'A'],
    ['G', 'R', 'Q', 'Y', 'K', 'V', 'Y', 'Y', 'H', 'Q'],
    ['Y', 'G', 'R', 'N', 'K', 'Y', 'R', 'H', 'Y', 'E'],
    ['J', 'E', 'M', 'A', 'S', 'L', 'F', 'I', 'Q', 'N']
]

result = find_largest_square(matrix)
print(result)
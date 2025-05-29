def find_largest_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    max_size = 0
    bottom_right_corner = (-1, -1)

    # Iterate over each possible starting point for a square
    for i in range(rows - 2):
        for j in range(cols - 2):
            # Check for squares of increasing size
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
    ['P', 'T', 'R', 'Z', 'B', 'V'],
    ['P', 'L', 'R', 'U', 'Q', 'U'],
    ['X', 'M', 'C', 'R', 'D', 'Y'],
    ['N', 'E', 'A', 'S', 'W', 'X'],
    ['C', 'I', 'L', 'L', 'L', 'Y'],
    ['J', 'F', 'L', 'L', 'L', 'Y'],
    ['H', 'G', 'L', 'L', 'L', 'E'],
    ['R', 'N', 'M', 'N', 'P', 'N'],
    ['W', 'F', 'L', 'M', 'L', 'E'],
    ['Z', 'W', 'C', 'Z', 'J', 'I']
]

result = find_largest_square(matrix)
print(result)
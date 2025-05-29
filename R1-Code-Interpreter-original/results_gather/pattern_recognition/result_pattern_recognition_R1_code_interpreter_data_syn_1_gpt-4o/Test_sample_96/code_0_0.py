def find_largest_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    max_size = 0
    bottom_right_corner = (-1, -1)

    # Start checking for squares of size 3x3 and larger
    for i in range(rows):
        for j in range(cols):
            # Check possible square sizes starting from 3
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
    ['N', 'G', 'U', 'O', 'U', 'K', 'G', 'C', 'Z', 'N', 'Z'],
    ['F', 'H', 'V', 'R', 'D', 'X', 'T', 'N', 'P', 'V', 'T'],
    ['B', 'B', 'B', 'Y', 'T', 'T', 'D', 'B', 'F', 'Z', 'F'],
    ['B', 'B', 'B', 'S', 'Q', 'U', 'K', 'M', 'R', 'H', 'C'],
    ['B', 'B', 'B', 'E', 'H', 'V', 'U', 'U', 'H', 'N', 'B'],
    ['F', 'O', 'S', 'R', 'Y', 'S', 'U', 'C', 'K', 'N', 'N'],
    ['I', 'X', 'J', 'G', 'S', 'Y', 'B', 'K', 'M', 'E', 'B']
]

result = find_largest_square(matrix)
print(result)
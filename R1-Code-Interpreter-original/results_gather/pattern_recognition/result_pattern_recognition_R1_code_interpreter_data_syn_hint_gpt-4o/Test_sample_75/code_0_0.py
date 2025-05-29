def find_largest_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    max_size = 0
    bottom_right_corner = (-1, -1)

    for i in range(rows):
        for j in range(cols):
            # Check for possible square sizes starting from 3
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
    ['J', 'K', 'S', 'P', 'X', 'L'],
    ['Q', 'V', 'G', 'W', 'W', 'W'],
    ['T', 'U', 'S', 'W', 'W', 'W'],
    ['N', 'X', 'N', 'W', 'W', 'W'],
    ['E', 'I', 'L', 'H', 'A', 'B'],
    ['R', 'Z', 'G', 'U', 'D', 'V'],
    ['J', 'S', 'B', 'F', 'U', 'X'],
    ['Y', 'K', 'O', 'V', 'U', 'E'],
    ['Z', 'F', 'X', 'R', 'I', 'N']
]

result = find_largest_square(matrix)
print(result)
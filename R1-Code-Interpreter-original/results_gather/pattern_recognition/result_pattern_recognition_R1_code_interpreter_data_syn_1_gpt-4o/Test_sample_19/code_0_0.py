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
    ['Q', 'L', 'H', 'E', 'X', 'Q', 'Y', 'D', 'Y', 'B', 'E', 'F'],
    ['L', 'I', 'A', 'M', 'B', 'G', 'X', 'C', 'N', 'J', 'V', 'W'],
    ['C', 'Q', 'W', 'N', 'F', 'A', 'B', 'U', 'O', 'P', 'B', 'X'],
    ['F', 'O', 'K', 'I', 'L', 'A', 'Z', 'F', 'L', 'L', 'L', 'G'],
    ['I', 'V', 'S', 'Q', 'N', 'T', 'M', 'C', 'L', 'L', 'L', 'E'],
    ['G', 'T', 'V', 'P', 'K', 'H', 'G', 'L', 'L', 'L', 'L', 'I']
]

result = find_largest_square(matrix)
print(result)
def find_largest_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    max_size = 0
    bottom_right_corner = None

    for i in range(rows):
        for j in range(cols):
            # Try to form a square with (i, j) as the top-left corner
            for size in range(3, min(rows - i, cols - j) + 1):
                # Check if all characters in the square are the same
                char = matrix[i][j]
                valid_square = True
                for x in range(i, i + size):
                    for y in range(j, j + size):
                        if matrix[x][y] != char:
                            valid_square = False
                            break
                    if not valid_square:
                        break
                
                if valid_square and size > max_size:
                    max_size = size
                    bottom_right_corner = (i + size - 1, j + size - 1)

    return bottom_right_corner

matrix = [
    ['S', 'P', 'P', 'P', 'P', 'W', 'V', 'L', 'Z', 'N', 'W', 'Y'],
    ['P', 'P', 'P', 'P', 'P', 'Q', 'G', 'C', 'F', 'C', 'S', 'U'],
    ['G', 'P', 'P', 'P', 'P', 'F', 'D', 'K', 'V', 'W', 'W', 'Y'],
    ['A', 'P', 'P', 'P', 'P', 'U', 'F', 'C', 'I', 'D', 'W', 'Q'],
    ['A', 'Z', 'E', 'B', 'Z', 'A', 'W', 'A', 'E', 'W', 'C', 'E'],
    ['C', 'V', 'L', 'X', 'R', 'G', 'W', 'H', 'J', 'N', 'T', 'J'],
    ['F', 'F', 'H', 'R', 'P', 'A', 'J', 'F', 'L', 'S', 'U', 'W'],
    ['Z', 'T', 'C', 'V', 'D', 'U', 'J', 'E', 'A', 'V', 'B', 'B'],
    ['A', 'L', 'T', 'Z', 'Q', 'F', 'S', 'J', 'A', 'J', 'G', 'J']
]

result = find_largest_square(matrix)
print(result)
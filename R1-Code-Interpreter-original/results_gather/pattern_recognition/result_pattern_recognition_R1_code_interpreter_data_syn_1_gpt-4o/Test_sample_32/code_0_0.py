def find_largest_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    max_size = 0
    bottom_right_corner = None

    for i in range(rows):
        for j in range(cols):
            # Try to form a square with top-left corner at (i, j)
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
    ['F', 'P', 'F', 'Y', 'F', 'H', 'X', 'U', 'E', 'W', 'R', 'Z', 'C'],
    ['O', 'M', 'G', 'M', 'J', 'E', 'A', 'G', 'R', 'B', 'C', 'Z', 'U'],
    ['K', 'Z', 'Q', 'Z', 'K', 'F', 'F', 'X', 'Q', 'S', 'B', 'B', 'B'],
    ['D', 'W', 'W', 'A', 'G', 'U', 'U', 'M', 'X', 'Q', 'B', 'B', 'B'],
    ['Y', 'W', 'C', 'F', 'Q', 'F', 'X', 'Q', 'W', 'X', 'B', 'B', 'B'],
    ['R', 'N', 'V', 'R', 'P', 'K', 'Q', 'Q', 'S', 'W', 'M', 'W', 'I'],
    ['G', 'Q', 'N', 'A', 'E', 'H', 'R', 'B', 'O', 'U', 'N', 'N', 'J']
]

result = find_largest_square(matrix)
print(result)
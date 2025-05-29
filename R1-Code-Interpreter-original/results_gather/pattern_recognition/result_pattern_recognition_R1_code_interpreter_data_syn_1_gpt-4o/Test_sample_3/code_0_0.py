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
    ['B', 'K', 'I', 'O', 'F', 'H', 'B', 'J', 'L'],
    ['F', 'J', 'L', 'G', 'P', 'D', 'S', 'Z', 'U'],
    ['I', 'L', 'Z', 'Z', 'G', 'W', 'Y', 'V', 'F'],
    ['J', 'T', 'D', 'L', 'Q', 'X', 'R', 'F', 'Z'],
    ['E', 'X', 'O', 'I', 'P', 'X', 'V', 'G', 'O'],
    ['K', 'F', 'R', 'Q', 'F', 'H', 'W', 'L', 'O'],
    ['S', 'H', 'W', 'C', 'F', 'Z', 'E', 'L', 'Z'],
    ['N', 'U', 'B', 'B', 'B', 'B', 'B', 'W', 'L'],
    ['S', 'B', 'B', 'B', 'B', 'B', 'B', 'N', 'T'],
    ['J', 'B', 'A', 'B', 'B', 'B', 'B', 'T', 'F'],
    ['K', 'X', 'O', 'B', 'B', 'B', 'B', 'O', 'B'],
    ['M', 'I', 'R', 'R', 'L', 'L', 'Q', 'Y', 'Q'],
    ['H', 'D', 'J', 'A', 'D', 'S', 'O', 'V', 'H']
]

result = find_largest_square(matrix)
print(result)
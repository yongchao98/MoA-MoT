def find_largest_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    max_size = 0
    bottom_right_corner = (-1, -1)

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
    ['M', 'F', 'Q', 'G', 'Q', 'B', 'D', 'Z', 'X', 'D', 'N', 'J'],
    ['W', 'W', 'I', 'D', 'M', 'O', 'B', 'U', 'S', 'S', 'S', 'Q'],
    ['F', 'N', 'V', 'I', 'F', 'S', 'T', 'P', 'S', 'S', 'S', 'V'],
    ['V', 'E', 'Y', 'B', 'Y', 'V', 'K', 'H', 'S', 'S', 'S', 'H'],
    ['I', 'Q', 'D', 'W', 'Q', 'L', 'D', 'R', 'P', 'B', 'D', 'O'],
    ['Y', 'C', 'G', 'D', 'E', 'G', 'O', 'N', 'I', 'D', 'D', 'R'],
    ['R', 'Q', 'K', 'R', 'C', 'I', 'W', 'G', 'S', 'P', 'A', 'H']
]

result = find_largest_square(matrix)
print(result)
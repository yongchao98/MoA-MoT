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
    ['D', 'D', 'L', 'C', 'F', 'G', 'U', 'H', 'H', 'P', 'Q', 'O'],
    ['P', 'N', 'F', 'A', 'A', 'A', 'C', 'H', 'A', 'Z', 'K', 'J'],
    ['V', 'K', 'L', 'A', 'A', 'A', 'K', 'U', 'K', 'M', 'L', 'O'],
    ['N', 'G', 'R', 'A', 'A', 'A', 'E', 'E', 'H', 'F', 'Z', 'W'],
    ['I', 'A', 'D', 'Z', 'C', 'E', 'G', 'D', 'E', 'A', 'G', 'A'],
    ['H', 'H', 'H', 'T', 'N', 'L', 'S', 'V', 'F', 'W', 'T', 'V'],
    ['L', 'T', 'U', 'J', 'T', 'B', 'W', 'I', 'Q', 'Q', 'N', 'S']
]

result = find_largest_square(matrix)
print(result)
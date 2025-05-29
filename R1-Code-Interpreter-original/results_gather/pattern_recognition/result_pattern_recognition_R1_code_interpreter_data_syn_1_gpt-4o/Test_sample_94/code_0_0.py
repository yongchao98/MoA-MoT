def find_largest_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    max_size = 0
    bottom_right_corner = None

    for i in range(rows):
        for j in range(cols):
            # Check for squares with top-left corner at (i, j)
            current_char = matrix[i][j]
            size = 2  # Start checking from size 3x3
            while i + size < rows and j + size < cols:
                # Check if all elements in the square are the same
                valid_square = True
                for x in range(i, i + size + 1):
                    for y in range(j, j + size + 1):
                        if matrix[x][y] != current_char:
                            valid_square = False
                            break
                    if not valid_square:
                        break
                
                if valid_square:
                    if size + 1 > max_size:
                        max_size = size + 1
                        bottom_right_corner = (i + size, j + size)
                else:
                    break
                size += 1

    return bottom_right_corner

matrix = [
    ['G', 'R', 'Z', 'K', 'I', 'V'],
    ['C', 'J', 'H', 'B', 'Y', 'K'],
    ['B', 'L', 'L', 'Q', 'U', 'X'],
    ['U', 'O', 'R', 'Q', 'T', 'R'],
    ['Z', 'E', 'L', 'O', 'O', 'O'],
    ['X', 'J', 'J', 'O', 'O', 'O'],
    ['C', 'Y', 'M', 'O', 'O', 'O']
]

result = find_largest_square(matrix)
print(result)
def find_largest_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    max_size = 0
    bottom_right_corner = None

    for i in range(rows):
        for j in range(cols):
            # Try to form a square with top-left corner at (i, j)
            current_char = matrix[i][j]
            size = 2  # Start checking from 3x3 squares
            while i + size < rows and j + size < cols:
                # Check if all characters in the square are the same
                valid_square = True
                for x in range(i, i + size + 1):
                    if not valid_square:
                        break
                    for y in range(j, j + size + 1):
                        if matrix[x][y] != current_char:
                            valid_square = False
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
    ['E', 'N', 'W', 'Y', 'F', 'A', 'R'],
    ['O', 'C', 'Z', 'U', 'J', 'D', 'Y'],
    ['Z', 'C', 'C', 'C', 'C', 'C', 'P'],
    ['N', 'J', 'S', 'C', 'C', 'C', 'A'],
    ['M', 'H', 'H', 'C', 'C', 'C', 'S'],
    ['P', 'W', 'A', 'A', 'P', 'M', 'U'],
    ['R', 'P', 'S', 'F', 'M', 'P', 'B'],
    ['H', 'T', 'B', 'C', 'X', 'J', 'B'],
    ['O', 'Z', 'X', 'A', 'P', 'E', 'V']
]

result = find_largest_square(matrix)
print(result)
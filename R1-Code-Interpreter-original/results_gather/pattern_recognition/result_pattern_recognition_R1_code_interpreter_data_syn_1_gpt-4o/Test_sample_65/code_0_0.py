def find_largest_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    max_size = 0
    bottom_right_corner = (-1, -1)

    for i in range(rows):
        for j in range(cols):
            # Try to form a square with (i, j) as the top-left corner
            current_char = matrix[i][j]
            size = 2  # Start checking from size 3x3
            while i + size < rows and j + size < cols:
                # Check if all characters in the square are the same
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
    ['E', 'Z', 'I', 'I', 'I', 'M', 'Y', 'V', 'F', 'O', 'X', 'U', 'Q', 'Y'],
    ['B', 'Y', 'I', 'I', 'I', 'Z', 'Y', 'F', 'X', 'G', 'G', 'Z', 'M', 'I'],
    ['L', 'C', 'I', 'I', 'I', 'J', 'L', 'M', 'S', 'U', 'G', 'E', 'T', 'X'],
    ['B', 'S', 'T', 'R', 'E', 'U', 'S', 'A', 'O', 'K', 'M', 'B', 'X', 'A'],
    ['N', 'F', 'T', 'S', 'C', 'A', 'E', 'C', 'W', 'Q', 'V', 'T', 'A', 'H'],
    ['T', 'D', 'V', 'U', 'H', 'M', 'E', 'H', 'S', 'Q', 'Q', 'R', 'Z', 'O']
]

result = find_largest_square(matrix)
print(result)
def find_largest_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    max_size = 0
    bottom_right_corner = None

    for i in range(rows):
        for j in range(cols):
            # Try to form a square with top-left corner at (i, j)
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
    ['C', 'S', 'A', 'F', 'O', 'U', 'I', 'T', 'C'],
    ['Y', 'Z', 'N', 'H', 'Y', 'B', 'B', 'G', 'O'],
    ['S', 'J', 'C', 'X', 'H', 'C', 'T', 'Q', 'I'],
    ['K', 'S', 'W', 'F', 'C', 'S', 'P', 'G', 'X'],
    ['G', 'G', 'G', 'G', 'Y', 'Q', 'A', 'Y', 'G'],
    ['G', 'G', 'G', 'G', 'M', 'S', 'M', 'K', 'M'],
    ['G', 'G', 'G', 'G', 'S', 'O', 'P', 'G', 'H'],
    ['G', 'G', 'G', 'G', 'G', 'A', 'Y', 'X', 'X'],
    ['D', 'X', 'V', 'U', 'X', 'R', 'A', 'L', 'T'],
    ['T', 'B', 'M', 'F', 'J', 'B', 'H', 'W', 'Q'],
    ['B', 'N', 'Z', 'K', 'Y', 'T', 'D', 'E', 'Q'],
    ['Y', 'E', 'T', 'Z', 'T', 'M', 'E', 'Z', 'K'],
    ['W', 'A', 'R', 'T', 'V', 'X', 'Z', 'H', 'X']
]

result = find_largest_square(matrix)
print(result)
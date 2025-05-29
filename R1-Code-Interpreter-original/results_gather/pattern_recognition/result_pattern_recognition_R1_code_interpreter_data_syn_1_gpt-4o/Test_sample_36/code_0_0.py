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
    ['S', 'P', 'W', 'I', 'W', 'H', 'N'],
    ['U', 'L', 'S', 'P', 'N', 'N', 'N'],
    ['M', 'V', 'H', 'H', 'N', 'N', 'N'],
    ['I', 'R', 'I', 'F', 'N', 'N', 'N'],
    ['W', 'M', 'I', 'D', 'F', 'P', 'S'],
    ['T', 'H', 'X', 'M', 'U', 'K', 'C'],
    ['Y', 'A', 'W', 'N', 'L', 'P', 'H'],
    ['O', 'V', 'I', 'S', 'N', 'Q', 'T'],
    ['T', 'K', 'O', 'S', 'E', 'L', 'I'],
    ['B', 'Y', 'Y', 'G', 'N', 'B', 'S'],
    ['U', 'Y', 'E', 'E', 'M', 'W', 'J'],
    ['W', 'O', 'C', 'C', 'V', 'L', 'W']
]

result = find_largest_square(matrix)
print(result)
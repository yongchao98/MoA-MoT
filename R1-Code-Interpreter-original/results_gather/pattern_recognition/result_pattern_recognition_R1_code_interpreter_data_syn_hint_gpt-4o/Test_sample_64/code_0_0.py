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
    ['F', 'N', 'U', 'V', 'P', 'Q'],
    ['C', 'I', 'I', 'I', 'V', 'D'],
    ['Y', 'I', 'I', 'I', 'M', 'R'],
    ['Z', 'I', 'I', 'I', 'S', 'J'],
    ['C', 'N', 'C', 'H', 'L', 'O'],
    ['O', 'S', 'A', 'Z', 'Y', 'H'],
    ['C', 'E', 'J', 'R', 'Y', 'Z'],
    ['X', 'A', 'V', 'N', 'P', 'E'],
    ['O', 'K', 'D', 'Q', 'Q', 'I'],
    ['E', 'H', 'G', 'N', 'P', 'F']
]

result = find_largest_square(matrix)
print(result)
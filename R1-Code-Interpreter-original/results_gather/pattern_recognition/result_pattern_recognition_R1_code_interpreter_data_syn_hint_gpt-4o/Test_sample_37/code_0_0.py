def find_largest_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    max_size = 0
    bottom_right_corner = None

    for i in range(rows):
        for j in range(cols):
            # Try to form squares of increasing size
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
    ['W', 'Y', 'D', 'J', 'V', 'K', 'C', 'E', 'B', 'H', 'Q', 'T'],
    ['M', 'B', 'O', 'C', 'M', 'P', 'P', 'A', 'U', 'W', 'P', 'L'],
    ['S', 'V', 'Y', 'Z', 'U', 'Z', 'U', 'K', 'W', 'R', 'S', 'Q'],
    ['O', 'Y', 'K', 'D', 'Y', 'C', 'A', 'W', 'V', 'N', 'U', 'R'],
    ['C', 'U', 'V', 'I', 'W', 'U', 'A', 'B', 'W', 'D', 'V', 'Y'],
    ['E', 'I', 'L', 'E', 'V', 'F', 'X', 'J', 'L', 'M', 'R', 'F'],
    ['D', 'K', 'D', 'D', 'D', 'D', 'D', 'N', 'R', 'L', 'X', 'H'],
    ['N', 'W', 'D', 'D', 'D', 'D', 'D', 'P', 'O', 'W', 'B', 'K'],
    ['D', 'F', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'A', 'E', 'W'],
    ['E', 'Z', 'D', 'D', 'D', 'D', 'D', 'H', 'I', 'I', 'S', 'D'],
    ['V', 'Q', 'D', 'D', 'D', 'D', 'D', 'M', 'R', 'O', 'K', 'P'],
    ['K', 'J', 'I', 'E', 'Q', 'L', 'F', 'T', 'S', 'I', 'A', 'D'],
    ['K', 'E', 'T', 'H', 'L', 'A', 'E', 'C', 'S', 'P', 'G', 'B'],
    ['H', 'L', 'Z', 'F', 'C', 'M', 'G', 'X', 'Z', 'L', 'O', 'F']
]

result = find_largest_square(matrix)
print(result)
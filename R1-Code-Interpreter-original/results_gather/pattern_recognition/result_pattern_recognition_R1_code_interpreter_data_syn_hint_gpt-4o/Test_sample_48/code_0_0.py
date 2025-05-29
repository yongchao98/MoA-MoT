def find_largest_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    max_size = 0
    bottom_right_corner = None

    for i in range(rows):
        for j in range(cols):
            # Check for squares of increasing size starting from (i, j)
            for size in range(3, min(rows - i, cols - j) + 1):
                # Check if all characters in the square are the same
                char = matrix[i][j]
                is_square = True
                for x in range(i, i + size):
                    for y in range(j, j + size):
                        if matrix[x][y] != char:
                            is_square = False
                            break
                    if not is_square:
                        break
                
                if is_square and size > max_size:
                    max_size = size
                    bottom_right_corner = (i + size - 1, j + size - 1)

    return bottom_right_corner

matrix = [
    ['J', 'U', 'O', 'R', 'Z', 'U', 'T', 'O'],
    ['M', 'U', 'Y', 'I', 'W', 'U', 'Y', 'C'],
    ['V', 'D', 'O', 'L', 'R', 'Q', 'S', 'E'],
    ['K', 'N', 'Z', 'L', 'O', 'S', 'T', 'E'],
    ['B', 'C', 'F', 'Z', 'D', 'L', 'C', 'W'],
    ['K', 'K', 'P', 'P', 'P', 'A', 'Z', 'T'],
    ['O', 'T', 'P', 'P', 'P', 'L', 'R', 'W'],
    ['F', 'Q', 'P', 'P', 'P', 'I', 'V', 'E'],
    ['Q', 'A', 'V', 'X', 'O', 'F', 'N', 'C'],
    ['S', 'L', 'K', 'B', 'U', 'Z', 'H', 'R'],
    ['A', 'W', 'L', 'R', 'M', 'S', 'N', 'N'],
    ['A', 'L', 'E', 'Q', 'L', 'V', 'P', 'O'],
    ['G', 'T', 'L', 'K', 'C', 'H', 'I', 'H'],
    ['U', 'T', 'Y', 'K', 'N', 'K', 'B', 'J']
]

result = find_largest_square(matrix)
print(result)
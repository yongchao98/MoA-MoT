def find_largest_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    max_size = 0
    bottom_right_corner = None

    for i in range(rows):
        for j in range(cols):
            # Check for squares of increasing size
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
    ['O', 'I', 'Y', 'C', 'Z', 'J', 'Y', 'V', 'O', 'F', 'R', 'X', 'T', 'A'],
    ['V', 'Y', 'V', 'I', 'F', 'M', 'H', 'R', 'Z', 'Q', 'Z', 'X', 'B', 'M'],
    ['I', 'B', 'R', 'R', 'O', 'Z', 'A', 'L', 'E', 'Z', 'D', 'B', 'G', 'X'],
    ['M', 'T', 'V', 'B', 'K', 'D', 'I', 'T', 'K', 'K', 'K', 'K', 'K', 'Y'],
    ['F', 'X', 'K', 'W', 'D', 'M', 'I', 'W', 'K', 'K', 'K', 'K', 'K', 'C'],
    ['Q', 'K', 'O', 'H', 'V', 'Q', 'P', 'T', 'K', 'K', 'K', 'K', 'K', 'A'],
    ['Q', 'I', 'B', 'V', 'L', 'H', 'K', 'A', 'K', 'K', 'K', 'K', 'K', 'M'],
    ['C', 'A', 'H', 'C', 'M', 'P', 'L', 'Z', 'K', 'K', 'K', 'K', 'K', 'D'],
    ['W', 'I', 'I', 'U', 'N', 'N', 'D', 'S', 'V', 'R', 'C', 'E', 'L', 'B'],
    ['P', 'B', 'S', 'U', 'N', 'Z', 'M', 'Y', 'V', 'U', 'Y', 'L', 'Q', 'A'],
    ['D', 'S', 'F', 'J', 'O', 'D', 'P', 'W', 'K', 'Q', 'P', 'J', 'V', 'M'],
    ['V', 'U', 'M', 'F', 'N', 'A', 'S', 'C', 'I', 'X', 'I', 'T', 'B', 'I']
]

result = find_largest_square(matrix)
print(result)
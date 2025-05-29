def find_largest_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    max_size = 0
    bottom_right_corner = None

    for i in range(rows):
        for j in range(cols):
            # Check for squares of increasing size starting from 3x3
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
    ['A', 'O', 'B', 'P', 'P', 'A', 'X', 'C', 'A', 'D', 'Q', 'J', 'F'],
    ['I', 'Z', 'T', 'H', 'D', 'V', 'T', 'T', 'Q', 'M', 'T', 'W', 'S'],
    ['U', 'L', 'H', 'O', 'E', 'R', 'K', 'P', 'Z', 'P', 'C', 'H', 'A'],
    ['S', 'H', 'P', 'N', 'Z', 'B', 'D', 'P', 'C', 'K', 'G', 'B', 'T'],
    ['G', 'F', 'O', 'O', 'O', 'Q', 'K', 'A', 'O', 'U', 'O', 'F', 'E'],
    ['P', 'C', 'O', 'O', 'O', 'G', 'M', 'H', 'Y', 'J', 'R', 'J', 'K'],
    ['G', 'L', 'O', 'O', 'O', 'I', 'I', 'J', 'W', 'W', 'P', 'I', 'Q']
]

result = find_largest_square(matrix)
print(result)
def find_largest_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    max_size = 0
    bottom_right_corner = None

    for i in range(rows):
        for j in range(cols):
            # Check for potential squares with side length >= 3
            for size in range(3, min(rows - i, cols - j) + 1):
                # Check if all elements in the square are the same
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
    ['T', 'Q', 'P', 'W', 'M', 'D', 'P', 'W'],
    ['O', 'S', 'C', 'T', 'A', 'O', 'Y', 'B'],
    ['F', 'M', 'Q', 'Q', 'L', 'I', 'S', 'R'],
    ['K', 'Y', 'D', 'P', 'B', 'W', 'E', 'A'],
    ['O', 'I', 'C', 'U', 'X', 'S', 'H', 'I'],
    ['R', 'H', 'H', 'H', 'C', 'L', 'T', 'I'],
    ['K', 'H', 'H', 'H', 'G', 'Z', 'E', 'W'],
    ['X', 'H', 'H', 'H', 'L', 'Z', 'M', 'D'],
    ['S', 'W', 'O', 'X', 'E', 'V', 'H', 'Y'],
    ['R', 'S', 'I', 'W', 'J', 'O', 'I', 'Z'],
    ['W', 'V', 'O', 'X', 'X', 'N', 'W', 'S'],
    ['L', 'D', 'X', 'Q', 'O', 'T', 'C', 'N']
]

result = find_largest_square(matrix)
print(result)
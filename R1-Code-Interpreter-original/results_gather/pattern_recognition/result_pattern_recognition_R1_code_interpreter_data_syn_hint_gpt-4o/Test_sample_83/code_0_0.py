def find_largest_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    max_side = 0
    bottom_right_corner = None

    for i in range(rows):
        for j in range(cols):
            # Check for potential squares with side length >= 3
            for side in range(3, min(rows - i, cols - j) + 1):
                # Check if all elements in the square are the same
                char = matrix[i][j]
                is_square = True
                for x in range(i, i + side):
                    for y in range(j, j + side):
                        if matrix[x][y] != char:
                            is_square = False
                            break
                    if not is_square:
                        break
                if is_square and side > max_side:
                    max_side = side
                    bottom_right_corner = (i + side - 1, j + side - 1)

    return bottom_right_corner

matrix = [
    ['T', 'W', 'V', 'B', 'X', 'J', 'K', 'X', 'O', 'X', 'O', 'L', 'Q'],
    ['I', 'T', 'T', 'I', 'S', 'N', 'S', 'T', 'N', 'P', 'L', 'H', 'H'],
    ['E', 'W', 'R', 'A', 'K', 'B', 'T', 'F', 'X', 'N', 'V', 'Z', 'X'],
    ['M', 'Y', 'R', 'R', 'R', 'H', 'Q', 'E', 'C', 'F', 'P', 'N', 'B'],
    ['R', 'A', 'R', 'R', 'R', 'K', 'I', 'C', 'Q', 'P', 'E', 'E', 'Q'],
    ['V', 'O', 'R', 'R', 'R', 'O', 'Z', 'O', 'E', 'P', 'H', 'X', 'X'],
    ['W', 'W', 'S', 'Q', 'W', 'A', 'K', 'P', 'Y', 'B', 'Z', 'R', 'U'],
    ['F', 'K', 'R', 'W', 'Q', 'E', 'N', 'O', 'W', 'K', 'F', 'Z', 'R'],
    ['B', 'Z', 'X', 'W', 'V', 'L', 'C', 'X', 'G', 'I', 'H', 'X', 'F']
]

result = find_largest_square(matrix)
print(result)
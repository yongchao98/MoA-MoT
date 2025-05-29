def find_largest_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    max_side = 0
    bottom_right_corner = None

    for i in range(rows - 2):
        for j in range(cols - 2):
            # Check for a 3x3 square
            if (matrix[i][j] == matrix[i][j+1] == matrix[i][j+2] ==
                matrix[i+1][j] == matrix[i+1][j+1] == matrix[i+1][j+2] ==
                matrix[i+2][j] == matrix[i+2][j+1] == matrix[i+2][j+2]):
                if 3 > max_side:
                    max_side = 3
                    bottom_right_corner = (i+2, j+2)

    return bottom_right_corner

matrix = [
    ['Q', 'S', 'I', 'K', 'M', 'K', 'S', 'W'],
    ['I', 'Q', 'B', 'C', 'M', 'P', 'M', 'P'],
    ['M', 'Q', 'W', 'J', 'O', 'C', 'C', 'C'],
    ['R', 'O', 'B', 'Q', 'E', 'C', 'C', 'C'],
    ['P', 'D', 'B', 'A', 'R', 'C', 'C', 'C'],
    ['Q', 'P', 'C', 'D', 'O', 'H', 'J', 'P']
]

result = find_largest_square(matrix)
print(result)
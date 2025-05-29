def find_largest_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    max_side = min(rows, cols)
    largest_square = (0, 0, 0)  # (side_length, bottom_right_row, bottom_right_col)

    for i in range(rows):
        for j in range(cols):
            char = matrix[i][j]
            for side in range(3, max_side + 1):
                if i + side > rows or j + side > cols:
                    break
                is_square = True
                for x in range(i, i + side):
                    for y in range(j, j + side):
                        if matrix[x][y] != char:
                            is_square = False
                            break
                    if not is_square:
                        break
                if is_square and side > largest_square[0]:
                    largest_square = (side, i + side - 1, j + side - 1)

    return [largest_square[1], largest_square[2]]

matrix = [
    ['G', 'B', 'J', 'O', 'T', 'E', 'J', 'K', 'Q', 'N', 'H', 'V'],
    ['I', 'G', 'D', 'S', 'U', 'I', 'H', 'P', 'N', 'V', 'J', 'G'],
    ['M', 'U', 'T', 'F', 'S', 'U', 'L', 'D', 'G', 'S', 'L', 'K'],
    ['S', 'B', 'F', 'N', 'R', 'L', 'F', 'Y', 'O', 'K', 'Q', 'M'],
    ['T', 'F', 'F', 'F', 'F', 'F', 'F', 'O', 'Q', 'Y', 'M', 'Z'],
    ['Z', 'F', 'F', 'F', 'F', 'F', 'F', 'O', 'W', 'U', 'E', 'C'],
    ['D', 'F', 'F', 'F', 'F', 'F', 'F', 'B', 'G', 'L', 'U', 'C'],
    ['M', 'F', 'F', 'F', 'F', 'F', 'F', 'X', 'Q', 'K', 'X', 'D'],
    ['O', 'F', 'F', 'F', 'F', 'F', 'F', 'P', 'A', 'M', 'A', 'U'],
    ['E', 'F', 'F', 'F', 'F', 'F', 'F', 'N', 'B', 'M', 'V', 'B'],
    ['N', 'N', 'U', 'X', 'U', 'G', 'U', 'S', 'S', 'N', 'Z', 'Y'],
    ['X', 'G', 'C', 'D', 'V', 'W', 'K', 'F', 'T', 'P', 'J', 'V'],
    ['M', 'F', 'K', 'Z', 'P', 'R', 'R', 'K', 'F', 'S', 'X', 'A']
]

result = find_largest_square(matrix)
print(result)
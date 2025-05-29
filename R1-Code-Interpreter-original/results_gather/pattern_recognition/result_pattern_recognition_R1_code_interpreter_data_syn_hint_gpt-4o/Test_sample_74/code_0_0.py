def find_largest_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    max_size = 0
    bottom_right_corner = (-1, -1)

    for i in range(rows):
        for j in range(cols):
            # Check for possible squares starting from size 3x3
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
    ['J', 'A', 'L', 'T', 'O', 'H'],
    ['Z', 'U', 'B', 'V', 'M', 'Z'],
    ['P', 'W', 'U', 'V', 'M', 'F'],
    ['S', 'R', 'Q', 'C', 'P', 'S'],
    ['N', 'Y', 'O', 'Q', 'T', 'M'],
    ['H', 'A', 'W', 'W', 'W', 'T'],
    ['K', 'S', 'W', 'W', 'W', 'P'],
    ['F', 'O', 'W', 'W', 'W', 'V'],
    ['L', 'Q', 'X', 'W', 'M', 'S'],
    ['U', 'O', 'O', 'W', 'Q', 'M'],
    ['I', 'I', 'Q', 'V', 'Y', 'M'],
    ['R', 'S', 'R', 'U', 'A', 'H']
]

result = find_largest_square(matrix)
print(result)
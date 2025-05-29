def find_largest_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    max_size = 0
    bottom_right_corner = (-1, -1)

    for i in range(rows):
        for j in range(cols):
            # Check for squares starting from size 3x3
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
    ['B', 'E', 'S', 'O', 'W', 'T', 'U', 'P', 'Q', 'J', 'F', 'F'],
    ['W', 'J', 'D', 'R', 'F', 'N', 'H', 'H', 'J', 'N', 'P', 'B'],
    ['S', 'D', 'O', 'N', 'N', 'N', 'N', 'O', 'R', 'A', 'T', 'X'],
    ['H', 'S', 'R', 'N', 'N', 'N', 'N', 'Y', 'E', 'Y', 'K', 'G'],
    ['C', 'L', 'L', 'N', 'N', 'N', 'N', 'O', 'Y', 'Y', 'P', 'H'],
    ['H', 'W', 'D', 'N', 'N', 'N', 'N', 'I', 'K', 'O', 'D', 'D'],
    ['P', 'P', 'V', 'R', 'R', 'H', 'Q', 'Q', 'F', 'A', 'N', 'E'],
    ['J', 'B', 'M', 'F', 'C', 'P', 'X', 'A', 'L', 'N', 'B', 'P'],
    ['Z', 'Y', 'O', 'A', 'Q', 'T', 'D', 'J', 'D', 'K', 'M', 'O']
]

result = find_largest_square(matrix)
print(result)
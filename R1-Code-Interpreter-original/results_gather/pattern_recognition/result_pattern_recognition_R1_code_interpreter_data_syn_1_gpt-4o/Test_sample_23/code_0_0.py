def find_largest_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    max_square_size = 0
    bottom_right_corner = None

    for i in range(rows):
        for j in range(cols):
            # Try to form a square with top-left corner at (i, j)
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
                
                if valid_square and size > max_square_size:
                    max_square_size = size
                    bottom_right_corner = (i + size - 1, j + size - 1)

    return bottom_right_corner

matrix = [
    ['N', 'Q', 'S', 'F', 'Z', 'W', 'T', 'N', 'M', 'S', 'H', 'G'],
    ['V', 'R', 'B', 'K', 'X', 'F', 'D', 'N', 'X', 'F', 'Q', 'S'],
    ['Z', 'B', 'N', 'U', 'L', 'I', 'V', 'C', 'V', 'G', 'W', 'B'],
    ['B', 'H', 'R', 'X', 'J', 'O', 'U', 'Q', 'E', 'O', 'T', 'G'],
    ['M', 'R', 'S', 'G', 'C', 'P', 'B', 'A', 'G', 'X', 'H', 'M'],
    ['A', 'R', 'K', 'N', 'G', 'G', 'P', 'A', 'E', 'R', 'Z', 'J'],
    ['J', 'A', 'K', 'Y', 'A', 'Z', 'I', 'V', 'E', 'R', 'B', 'X'],
    ['B', 'L', 'H', 'P', 'D', 'V', 'Y', 'X', 'N', 'B', 'G', 'K'],
    ['Q', 'U', 'B', 'G', 'S', 'Q', 'W', 'M', 'V', 'U', 'R', 'V'],
    ['B', 'C', 'L', 'L', 'H', 'X', 'K', 'B', 'U', 'H', 'T', 'N'],
    ['N', 'P', 'N', 'E', 'J', 'N', 'X', 'Z', 'J', 'S', 'L', 'K'],
    ['T', 'Q', 'J', 'D', 'T', 'W', 'N', 'K', 'Q', 'F', 'E', 'X'],
    ['K', 'X', 'I', 'C', 'O', 'L', 'S', 'M', 'M', 'M', 'B', 'Y'],
    ['X', 'J', 'W', 'Q', 'K', 'L', 'G', 'M', 'M', 'M', 'C', 'M'],
    ['V', 'D', 'T', 'N', 'X', 'L', 'I', 'M', 'M', 'M', 'T', 'A']
]

result = find_largest_square(matrix)
print(result)
def find_largest_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    max_size = 0
    bottom_right_corner = None

    for i in range(rows):
        for j in range(cols):
            # Try to form squares of increasing size starting from 3x3
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
                
                if valid_square and size > max_size:
                    max_size = size
                    bottom_right_corner = (i + size - 1, j + size - 1)

    return bottom_right_corner

matrix = [
    ['J', 'R', 'Z', 'T', 'F', 'B', 'R'],
    ['N', 'F', 'J', 'S', 'Q', 'U', 'D'],
    ['F', 'B', 'H', 'U', 'X', 'L', 'W'],
    ['Y', 'M', 'K', 'U', 'L', 'H', 'N'],
    ['C', 'U', 'H', 'I', 'J', 'B', 'Q'],
    ['I', 'N', 'N', 'N', 'T', 'Z', 'I'],
    ['A', 'N', 'N', 'N', 'E', 'D', 'H'],
    ['B', 'N', 'N', 'N', 'X', 'U', 'R'],
    ['O', 'X', 'O', 'P', 'Q', 'R', 'I']
]

result = find_largest_square(matrix)
print(result)
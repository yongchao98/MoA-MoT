def find_largest_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    max_size = 0
    bottom_right_corner = None

    for i in range(rows):
        for j in range(cols):
            # Try to form squares of increasing size starting from 3x3
            for size in range(3, min(rows - i, cols - j) + 1):
                # Check if the square of this size is valid
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
    ['T', 'N', 'G', 'F', 'I', 'M', 'O'],
    ['I', 'U', 'G', 'I', 'C', 'S', 'V'],
    ['I', 'N', 'U', 'Z', 'R', 'Y', 'C'],
    ['V', 'V', 'Y', 'A', 'G', 'O', 'R'],
    ['Y', 'N', 'Z', 'Y', 'D', 'M', 'V'],
    ['P', 'W', 'S', 'X', 'C', 'E', 'K'],
    ['Q', 'B', 'W', 'P', 'O', 'A', 'D'],
    ['N', 'Z', 'C', 'P', 'T', 'C', 'N'],
    ['J', 'A', 'M', 'Z', 'H', 'E', 'J'],
    ['S', 'Q', 'O', 'Z', 'B', 'P', 'B'],
    ['X', 'X', 'X', 'U', 'C', 'C', 'O'],
    ['X', 'X', 'X', 'O', 'J', 'U', 'R'],
    ['X', 'X', 'X', 'O', 'L', 'R', 'H']
]

result = find_largest_square(matrix)
print(result)
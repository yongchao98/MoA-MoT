def find_largest_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    max_size = 0
    bottom_right_corner = (-1, -1)

    # Iterate over each possible top-left corner of a square
    for i in range(rows):
        for j in range(cols):
            # Try to form squares of increasing size
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
                
                # If a valid square is found, update the largest square
                if is_square and size > max_size:
                    max_size = size
                    bottom_right_corner = (i + size - 1, j + size - 1)

    return bottom_right_corner

matrix = [
    ['Q', 'Q', 'N', 'P', 'I', 'O'],
    ['N', 'M', 'P', 'W', 'M', 'K'],
    ['F', 'A', 'W', 'M', 'W', 'S'],
    ['C', 'E', 'E', 'E', 'D', 'V'],
    ['J', 'E', 'E', 'E', 'F', 'Z'],
    ['W', 'E', 'E', 'E', 'N', 'I'],
    ['I', 'Y', 'K', 'L', 'X', 'J'],
    ['P', 'W', 'L', 'G', 'H', 'P'],
    ['Q', 'U', 'M', 'G', 'D', 'Q'],
    ['U', 'M', 'W', 'H', 'D', 'U'],
    ['F', 'U', 'M', 'X', 'A', 'S']
]

result = find_largest_square(matrix)
print(result)
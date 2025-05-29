def find_largest_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    max_size = 0
    bottom_right_corner = None

    for i in range(rows):
        for j in range(cols):
            # Try to form a square with (i, j) as the top-left corner
            current_char = matrix[i][j]
            for size in range(3, min(rows - i, cols - j) + 1):
                # Check if all characters in the square are the same
                valid_square = True
                for x in range(size):
                    for y in range(size):
                        if matrix[i + x][j + y] != current_char:
                            valid_square = False
                            break
                    if not valid_square:
                        break
                
                if valid_square and size > max_size:
                    max_size = size
                    bottom_right_corner = (i + size - 1, j + size - 1)

    return bottom_right_corner

matrix = [
    ['Q', 'W', 'Z', 'P', 'S', 'V', 'Q', 'C'],
    ['E', 'S', 'E', 'G', 'G', 'G', 'G', 'I'],
    ['A', 'S', 'M', 'G', 'G', 'G', 'G', 'J'],
    ['R', 'G', 'B', 'G', 'G', 'G', 'G', 'U'],
    ['A', 'A', 'X', 'G', 'G', 'G', 'G', 'M'],
    ['R', 'C', 'U', 'C', 'V', 'B', 'A', 'D'],
    ['P', 'T', 'R', 'Q', 'D', 'O', 'N', 'W'],
    ['R', 'B', 'E', 'X', 'O', 'L', 'M', 'M'],
    ['M', 'G', 'F', 'J', 'W', 'N', 'Q', 'N'],
    ['X', 'N', 'P', 'M', 'Y', 'N', 'X', 'D'],
    ['U', 'V', 'C', 'S', 'C', 'H', 'C', 'K'],
    ['L', 'X', 'C', 'W', 'U', 'J', 'K', 'F'],
    ['Q', 'J', 'M', 'V', 'Q', 'H', 'X', 'O'],
    ['K', 'O', 'C', 'E', 'S', 'R', 'P', 'E']
]

result = find_largest_square(matrix)
print(result)
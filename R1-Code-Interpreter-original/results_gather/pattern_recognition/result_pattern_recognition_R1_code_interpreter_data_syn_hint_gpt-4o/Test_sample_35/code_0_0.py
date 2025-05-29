def find_largest_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    max_size = 0
    bottom_right_corner = None

    for i in range(rows):
        for j in range(cols):
            # Check for possible square sizes starting from 3x3
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
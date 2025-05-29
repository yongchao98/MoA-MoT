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
                # Check if the square of this size is uniform
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
    ['C', 'M', 'O', 'G', 'B', 'C', 'W', 'G', 'L', 'P'],
    ['E', 'P', 'C', 'S', 'M', 'W', 'V', 'D', 'K', 'J'],
    ['L', 'O', 'L', 'V', 'Y', 'U', 'J', 'B', 'V', 'H'],
    ['Z', 'Q', 'K', 'T', 'I', 'A', 'Y', 'E', 'G', 'R'],
    ['H', 'H', 'Q', 'O', 'S', 'W', 'J', 'N', 'C', 'M'],
    ['Q', 'V', 'I', 'T', 'C', 'G', 'Y', 'S', 'X', 'O'],
    ['P', 'B', 'B', 'W', 'D', 'E', 'W', 'M', 'B', 'S'],
    ['A', 'L', 'R', 'M', 'P', 'C', 'C', 'M', 'A', 'Y'],
    ['P', 'C', 'Z', 'B', 'E', 'B', 'S', 'G', 'I', 'X'],
    ['B', 'F', 'Y', 'O', 'F', 'F', 'F', 'F', 'C', 'D'],
    ['V', 'K', 'E', 'G', 'F', 'F', 'F', 'F', 'Y', 'C'],
    ['G', 'D', 'R', 'P', 'F', 'F', 'F', 'F', 'Z', 'N'],
    ['V', 'M', 'F', 'M', 'F', 'F', 'F', 'F', 'T', 'W'],
    ['R', 'B', 'P', 'K', 'B', 'Y', 'M', 'N', 'L', 'Y'],
    ['Y', 'J', 'Z', 'Z', 'Q', 'P', 'X', 'R', 'P', 'Z']
]

result = find_largest_square(matrix)
print(result)
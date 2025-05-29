def find_largest_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    max_size = 0
    bottom_right_corner = (-1, -1)

    for i in range(rows):
        for j in range(cols):
            # Try to form a square with (i, j) as the top-left corner
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
    ['H', 'U', 'A', 'A', 'I', 'W', 'X', 'X', 'F', 'A', 'K', 'C', 'J', 'I'],
    ['Q', 'S', 'W', 'F', 'P', 'A', 'B', 'I', 'P', 'V', 'V', 'R', 'J', 'I'],
    ['T', 'N', 'Z', 'Z', 'U', 'X', 'M', 'C', 'B', 'G', 'B', 'G', 'X', 'U'],
    ['S', 'B', 'F', 'J', 'V', 'M', 'U', 'N', 'V', 'Y', 'L', 'B', 'N', 'P'],
    ['X', 'G', 'S', 'A', 'G', 'Q', 'V', 'K', 'B', 'E', 'W', 'L', 'J', 'B'],
    ['S', 'N', 'K', 'M', 'X', 'F', 'P', 'H', 'F', 'V', 'N', 'U', 'S', 'R'],
    ['I', 'O', 'R', 'N', 'R', 'W', 'F', 'K', 'R', 'P', 'H', 'K', 'N', 'D'],
    ['Y', 'Y', 'E', 'X', 'B', 'G', 'V', 'G', 'G', 'G', 'G', 'G', 'Q', 'Y'],
    ['M', 'W', 'G', 'E', 'N', 'T', 'F', 'G', 'G', 'G', 'G', 'G', 'B', 'Q'],
    ['S', 'G', 'D', 'S', 'Q', 'G', 'I', 'G', 'G', 'G', 'G', 'G', 'A', 'B'],
    ['Z', 'I', 'P', 'L', 'V', 'S', 'X', 'G', 'G', 'G', 'G', 'G', 'C', 'Z'],
    ['K', 'H', 'B', 'A', 'X', 'N', 'Y', 'G', 'G', 'G', 'G', 'G', 'T', 'O'],
    ['I', 'H', 'U', 'I', 'V', 'W', 'X', 'A', 'I', 'A', 'Y', 'Z', 'W', 'F'],
    ['G', 'B', 'Z', 'R', 'V', 'B', 'L', 'R', 'Q', 'X', 'P', 'G', 'U', 'O'],
    ['D', 'M', 'Z', 'G', 'Z', 'F', 'F', 'D', 'M', 'E', 'W', 'X', 'R', 'Q']
]

result = find_largest_square(matrix)
print(result)
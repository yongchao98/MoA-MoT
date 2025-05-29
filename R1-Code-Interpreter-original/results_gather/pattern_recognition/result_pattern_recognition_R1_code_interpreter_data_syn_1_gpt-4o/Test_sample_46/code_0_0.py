def find_largest_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    max_size = 0
    bottom_right_corner = (-1, -1)

    for i in range(rows):
        for j in range(cols):
            # Check for squares of increasing size starting from 3
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
    ['K', 'R', 'J', 'Z', 'M', 'P', 'Y', 'Y', 'Z', 'F', 'Y', 'W', 'Z', 'R', 'Y'],
    ['G', 'W', 'U', 'X', 'R', 'L', 'C', 'E', 'S', 'A', 'I', 'P', 'Q', 'I', 'I'],
    ['M', 'Q', 'J', 'C', 'Y', 'V', 'P', 'A', 'Y', 'R', 'F', 'E', 'S', 'T', 'N'],
    ['O', 'Y', 'D', 'A', 'B', 'O', 'B', 'E', 'B', 'A', 'X', 'W', 'Y', 'Y', 'Y'],
    ['K', 'N', 'N', 'M', 'Z', 'Q', 'F', 'X', 'F', 'E', 'D', 'F', 'Y', 'Y', 'Y'],
    ['A', 'N', 'K', 'Y', 'E', 'A', 'L', 'N', 'R', 'D', 'F', 'R', 'Y', 'Y', 'Y'],
    ['W', 'Z', 'P', 'W', 'P', 'C', 'C', 'X', 'C', 'F', 'S', 'T', 'X', 'Y', 'M'],
    ['N', 'B', 'M', 'E', 'A', 'M', 'D', 'E', 'K', 'H', 'K', 'V', 'V', 'C', 'J'],
    ['G', 'J', 'I', 'W', 'I', 'N', 'T', 'Q', 'O', 'H', 'C', 'F', 'M', 'B', 'S']
]

result = find_largest_square(matrix)
print(result)
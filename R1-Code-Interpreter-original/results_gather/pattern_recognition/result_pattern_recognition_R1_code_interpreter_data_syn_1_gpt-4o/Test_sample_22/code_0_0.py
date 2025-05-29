def find_largest_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    max_size = 0
    bottom_right_corner = (-1, -1)

    for i in range(rows):
        for j in range(cols):
            # Check for potential squares starting from (i, j)
            current_char = matrix[i][j]
            size = 1
            while (i + size < rows) and (j + size < cols):
                # Check if the square of size+1 is valid
                valid_square = True
                for k in range(size + 1):
                    if matrix[i + size][j + k] != current_char or matrix[i + k][j + size] != current_char:
                        valid_square = False
                        break
                if not valid_square:
                    break
                size += 1

            # Update the largest square found
            if size > max_size and size >= 3:
                max_size = size
                bottom_right_corner = (i + size - 1, j + size - 1)

    return bottom_right_corner

matrix = [
    ['K', 'X', 'I', 'R', 'T', 'M', 'J', 'C', 'W', 'A', 'Z'],
    ['G', 'U', 'P', 'M', 'T', 'U', 'E', 'V', 'G', 'I', 'M'],
    ['E', 'C', 'L', 'V', 'J', 'J', 'J', 'F', 'T', 'K', 'L'],
    ['A', 'J', 'M', 'W', 'J', 'J', 'J', 'Q', 'T', 'O', 'P'],
    ['C', 'E', 'U', 'X', 'J', 'J', 'J', 'X', 'A', 'U', 'G'],
    ['O', 'L', 'H', 'F', 'F', 'O', 'Y', 'Y', 'Z', 'O', 'K'],
    ['C', 'K', 'C', 'T', 'Z', 'G', 'E', 'G', 'M', 'L', 'D'],
    ['P', 'G', 'H', 'Y', 'U', 'C', 'U', 'Y', 'S', 'E', 'Y'],
    ['R', 'U', 'Z', 'H', 'R', 'C', 'E', 'Q', 'R', 'J', 'V'],
    ['R', 'Q', 'Y', 'T', 'K', 'X', 'I', 'V', 'U', 'I', 'S'],
    ['D', 'A', 'R', 'H', 'Q', 'N', 'H', 'A', 'B', 'I', 'Y'],
    ['T', 'K', 'X', 'F', 'D', 'D', 'G', 'A', 'T', 'M', 'I']
]

result = find_largest_square(matrix)
print(result)
def find_largest_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    max_size = min(rows, cols)
    
    for size in range(3, max_size + 1):
        for i in range(rows - size + 1):
            for j in range(cols - size + 1):
                char = matrix[i][j]
                is_square = True
                for k in range(size):
                    if not (all(matrix[i + k][j + l] == char for l in range(size)) and
                            all(matrix[i + l][j + k] == char for l in range(size))):
                        is_square = False
                        break
                if is_square:
                    bottom_right = (i + size - 1, j + size - 1)
    
    return bottom_right

matrix = [
    ['Q', 'I', 'A', 'Y', 'W', 'B', 'D', 'K', 'K', 'E', 'J', 'F'],
    ['Q', 'R', 'M', 'U', 'E', 'Q', 'I', 'I', 'K', 'S', 'X', 'K'],
    ['Y', 'D', 'R', 'E', 'F', 'U', 'I', 'W', 'H', 'V', 'Y', 'I'],
    ['L', 'Z', 'S', 'E', 'F', 'P', 'A', 'M', 'D', 'M', 'Y', 'U'],
    ['X', 'P', 'C', 'C', 'C', 'C', 'H', 'Z', 'V', 'O', 'V', 'V'],
    ['L', 'T', 'C', 'C', 'C', 'C', 'P', 'Y', 'P', 'A', 'Q', 'V'],
    ['X', 'O', 'C', 'C', 'C', 'C', 'P', 'C', 'S', 'X', 'O', 'B'],
    ['P', 'N', 'C', 'C', 'C', 'C', 'Z', 'E', 'A', 'X', 'B', 'L'],
    ['K', 'R', 'U', 'A', 'O', 'Z', 'C', 'I', 'G', 'S', 'L', 'E'],
    ['D', 'G', 'O', 'Q', 'O', 'S', 'I', 'U', 'P', 'Z', 'E', 'N']
]

bottom_right_corner = find_largest_square(matrix)
print(bottom_right_corner)
def find_largest_square(matrix):
    n = len(matrix)
    max_size = 0
    bottom_right_corner = (-1, -1)
    
    for i in range(n):
        for j in range(n):
            # Check for squares of increasing size
            for size in range(3, n - max(i, j) + 1):
                if all(matrix[i][j] == matrix[i + k][j] == matrix[i][j + k] == matrix[i + size - 1][j + k] == matrix[i + k][j + size - 1] for k in range(size)):
                    if size > max_size:
                        max_size = size
                        bottom_right_corner = (i + size - 1, j + size - 1)
    
    return bottom_right_corner

matrix = [
    ['L', 'D', 'H', 'X', 'Q', 'G', 'D', 'V', 'E', 'S'],
    ['E', 'H', 'Y', 'T', 'D', 'Y', 'R', 'T', 'N', 'M'],
    ['K', 'K', 'K', 'K', 'A', 'M', 'E', 'M', 'M', 'K'],
    ['K', 'K', 'K', 'K', 'E', 'W', 'X', 'D', 'L', 'N'],
    ['K', 'K', 'K', 'K', 'Z', 'D', 'U', 'D', 'T', 'N'],
    ['K', 'K', 'K', 'K', 'Y', 'N', 'V', 'J', 'V', 'F'],
    ['Z', 'C', 'T', 'J', 'M', 'G', 'V', 'D', 'X', 'P'],
    ['G', 'D', 'A', 'J', 'R', 'U', 'Y', 'H', 'K', 'L'],
    ['S', 'U', 'N', 'E', 'R', 'A', 'Y', 'X', 'N', 'Y'],
    ['K', 'Q', 'I', 'O', 'P', 'N', 'R', 'G', 'N', 'V']
]

result = find_largest_square(matrix)
print(result)
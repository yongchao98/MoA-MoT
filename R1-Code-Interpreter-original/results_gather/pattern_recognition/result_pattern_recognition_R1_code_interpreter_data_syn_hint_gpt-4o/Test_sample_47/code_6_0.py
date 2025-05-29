def find_largest_square(matrix):
    n = len(matrix)
    max_size = 0
    bottom_right_corner = (-1, -1)
    
    for i in range(n):
        for j in range(n):
            # Try to form squares of increasing size
            for l in range(3, n - max(i, j) + 1):
                # Check if all characters in the square are the same
                char = matrix[i][j]
                valid_square = True
                for x in range(i, i + l):
                    for y in range(j, j + l):
                        if matrix[x][y] != char:
                            valid_square = False
                            break
                    if not valid_square:
                        break
                
                if valid_square and l > max_size:
                    max_size = l
                    bottom_right_corner = (i + l - 1, j + l - 1)
    
    return bottom_right_corner

matrix = [
    ['O', 'C', 'D', 'G', 'G', 'Q', 'J', 'Q', 'U', 'K'],
    ['X', 'U', 'E', 'C', 'B', 'R', 'V', 'S', 'K', 'M'],
    ['E', 'G', 'Z', 'R', 'R', 'J', 'R', 'Y', 'O', 'T'],
    ['W', 'U', 'J', 'J', 'V', 'P', 'G', 'C', 'Z', 'X'],
    ['K', 'T', 'F', 'L', 'P', 'P', 'N', 'H', 'N', 'I'],
    ['K', 'K', 'K', 'Z', 'E', 'Q', 'E', 'Q', 'C', 'Y'],
    ['K', 'K', 'K', 'U', 'U', 'W', 'U', 'K', 'S', 'B'],
    ['K', 'K', 'K', 'O', 'B', 'O', 'B', 'F', 'S', 'V'],
    ['B', 'O', 'L', 'B', 'L', 'Z', 'T', 'Z', 'F', 'T'],
    ['N', 'B', 'T', 'X', 'O', 'Z', 'N', 'L', 'L', 'Y']
]

result = find_largest_square(matrix)
print(result)
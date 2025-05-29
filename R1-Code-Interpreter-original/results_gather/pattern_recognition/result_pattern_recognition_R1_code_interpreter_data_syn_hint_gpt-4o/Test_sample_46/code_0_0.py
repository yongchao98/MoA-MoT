def find_squares(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    max_square_size = min(rows, cols)
    
    for size in range(3, max_square_size + 1):
        for i in range(rows - size + 1):
            for j in range(cols - size + 1):
                char = matrix[i][j]
                is_square = True
                for x in range(size):
                    for y in range(size):
                        if matrix[i + x][j + y] != char:
                            is_square = False
                            break
                    if not is_square:
                        break
                if is_square:
                    return [i + size - 1, j + size - 1]

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

result = find_squares(matrix)
print(result)
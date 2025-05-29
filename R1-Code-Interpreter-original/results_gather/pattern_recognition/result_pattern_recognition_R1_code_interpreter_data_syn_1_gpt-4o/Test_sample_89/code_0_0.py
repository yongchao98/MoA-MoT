def find_largest_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    max_size = min(rows, cols)
    
    for size in range(3, max_size + 1):
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
    ['T', 'Q', 'P', 'W', 'M', 'D', 'P', 'W'],
    ['O', 'S', 'C', 'T', 'A', 'O', 'Y', 'B'],
    ['F', 'M', 'Q', 'Q', 'L', 'I', 'S', 'R'],
    ['K', 'Y', 'D', 'P', 'B', 'W', 'E', 'A'],
    ['O', 'I', 'C', 'U', 'X', 'S', 'H', 'I'],
    ['R', 'H', 'H', 'H', 'C', 'L', 'T', 'I'],
    ['K', 'H', 'H', 'H', 'G', 'Z', 'E', 'W'],
    ['X', 'H', 'H', 'H', 'L', 'Z', 'M', 'D'],
    ['S', 'W', 'O', 'X', 'E', 'V', 'H', 'Y'],
    ['R', 'S', 'I', 'W', 'J', 'O', 'I', 'Z'],
    ['W', 'V', 'O', 'X', 'X', 'N', 'W', 'S'],
    ['L', 'D', 'X', 'Q', 'O', 'T', 'C', 'N']
]

print(find_largest_square(matrix))
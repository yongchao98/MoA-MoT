def find_square(matrix):
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
    ['I', 'I', 'Z', 'Q', 'A', 'R', 'L', 'S', 'E', 'T', 'F', 'F', 'I', 'J', 'Z'],
    ['S', 'S', 'S', 'S', 'S', 'H', 'P', 'V', 'T', 'L', 'D', 'C', 'C', 'X', 'J'],
    ['E', 'S', 'S', 'S', 'S', 'Z', 'H', 'Q', 'G', 'H', 'Q', 'I', 'O', 'E', 'A'],
    ['R', 'S', 'S', 'S', 'S', 'E', 'Y', 'X', 'F', 'M', 'G', 'F', 'K', 'R', 'V'],
    ['X', 'S', 'S', 'S', 'S', 'A', 'K', 'O', 'M', 'Z', 'B', 'E', 'P', 'T', 'B'],
    ['A', 'Y', 'R', 'E', 'A', 'W', 'H', 'K', 'Z', 'Y', 'Q', 'E', 'Z', 'W', 'C'],
    ['N', 'F', 'B', 'E', 'C', 'K', 'H', 'M', 'W', 'A', 'L', 'R', 'X', 'P', 'K'],
    ['H', 'R', 'L', 'U', 'Q', 'B', 'G', 'I', 'U', 'E', 'Y', 'T', 'M', 'S', 'Z'],
    ['L', 'E', 'V', 'S', 'U', 'X', 'L', 'H', 'V', 'U', 'R', 'V', 'Q', 'O', 'O'],
    ['Q', 'I', 'T', 'C', 'T', 'C', 'S', 'V', 'B', 'O', 'B', 'I', 'F', 'B', 'B']
]

result = find_square(matrix)
print(result)
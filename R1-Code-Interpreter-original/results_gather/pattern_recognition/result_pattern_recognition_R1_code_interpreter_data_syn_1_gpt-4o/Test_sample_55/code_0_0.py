def find_squares(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    min_size = 3
    result = None

    for i in range(rows):
        for j in range(cols):
            # Check for squares of increasing size
            for size in range(min_size, min(rows - i, cols - j) + 1):
                # Check if all characters in the square are the same
                char = matrix[i][j]
                is_square = True
                for x in range(i, i + size):
                    for y in range(j, j + size):
                        if matrix[x][y] != char:
                            is_square = False
                            break
                    if not is_square:
                        break
                if is_square:
                    result = (i + size - 1, j + size - 1)
    
    return result

matrix = [
    ['K', 'A', 'L', 'I', 'T', 'Y', 'V', 'A', 'L', 'Q', 'F', 'I', 'Y', 'O', 'T'],
    ['X', 'I', 'Q', 'W', 'E', 'I', 'A', 'K', 'T', 'J', 'D', 'M', 'E', 'F', 'S'],
    ['S', 'I', 'Y', 'A', 'W', 'J', 'Y', 'V', 'V', 'V', 'J', 'R', 'J', 'L', 'R'],
    ['U', 'C', 'O', 'X', 'Q', 'T', 'M', 'V', 'V', 'V', 'X', 'M', 'E', 'Z', 'M'],
    ['T', 'H', 'E', 'R', 'Q', 'U', 'X', 'V', 'V', 'V', 'L', 'S', 'O', 'T', 'H'],
    ['B', 'R', 'O', 'T', 'X', 'J', 'J', 'I', 'W', 'D', 'J', 'X', 'U', 'S', 'D'],
    ['Q', 'U', 'D', 'Z', 'A', 'F', 'W', 'Z', 'F', 'E', 'K', 'T', 'L', 'S', 'N'],
    ['O', 'O', 'C', 'H', 'C', 'Z', 'M', 'H', 'R', 'D', 'C', 'P', 'T', 'I', 'Z'],
    ['B', 'W', 'S', 'C', 'D', 'U', 'T', 'D', 'Q', 'A', 'E', 'X', 'A', 'A', 'Y'],
    ['N', 'S', 'X', 'A', 'V', 'H', 'Z', 'H', 'L', 'W', 'R', 'P', 'V', 'E', 'R'],
    ['K', 'O', 'O', 'B', 'X', 'M', 'Y', 'D', 'H', 'G', 'D', 'E', 'V', 'T', 'S'],
    ['C', 'N', 'C', 'M', 'N', 'U', 'R', 'S', 'Y', 'T', 'F', 'M', 'Q', 'Q', 'F']
]

result = find_squares(matrix)
print(result)
def find_squares(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    min_size = 3
    result = []

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
                    result.append((i + size - 1, j + size - 1))
    
    return result

matrix = [
    ['I', 'E', 'L', 'F', 'L', 'B'],
    ['Q', 'V', 'L', 'Q', 'F', 'X'],
    ['T', 'P', 'K', 'J', 'U', 'R'],
    ['D', 'Q', 'D', 'W', 'A', 'Y'],
    ['W', 'S', 'S', 'S', 'S', 'S'],
    ['A', 'S', 'S', 'S', 'Y', 'E'],
    ['B', 'S', 'S', 'S', 'G', 'J'],
    ['V', 'R', 'T', 'I', 'X', 'M']
]

squares = find_squares(matrix)
print(squares)
def find_squares(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    min_size = 3
    result = []

    for i in range(rows):
        for j in range(cols):
            # Check for squares starting from size 3x3
            for size in range(min_size, min(rows - i, cols - j) + 1):
                char = matrix[i][j]
                is_square = True
                # Check if all characters in the square are the same
                for x in range(i, i + size):
                    for y in range(j, j + size):
                        if matrix[x][y] != char:
                            is_square = False
                            break
                    if not is_square:
                        break
                if is_square:
                    result.append((i + size - 1, j + size - 1))
    
    # Return the last found square's bottom-right corner
    return result[-1] if result else None

matrix = [
    ['P', 'S', 'K', 'Y', 'A', 'O', 'X', 'B', 'H', 'C'],
    ['T', 'Y', 'Z', 'V', 'T', 'M', 'K', 'X', 'A', 'T'],
    ['D', 'P', 'R', 'O', 'Z', 'K', 'K', 'K', 'Y', 'J'],
    ['R', 'X', 'K', 'W', 'M', 'K', 'K', 'K', 'N', 'G'],
    ['B', 'V', 'B', 'T', 'D', 'K', 'K', 'K', 'K', 'I'],
    ['B', 'M', 'I', 'D', 'X', 'M', 'P', 'E', 'A', 'K']
]

result = find_squares(matrix)
print(result)
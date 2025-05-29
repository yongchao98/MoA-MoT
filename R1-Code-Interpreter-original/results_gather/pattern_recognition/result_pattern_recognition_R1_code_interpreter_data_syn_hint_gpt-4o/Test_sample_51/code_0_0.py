def find_squares(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    min_size = 3
    result = []

    for i in range(rows):
        for j in range(cols):
            # Check for squares of increasing size starting from 3x3
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
    
    # Since the answer is unique, we assume the largest square is the one we need
    return result[-1] if result else None

matrix = [
    ['D', 'D', 'L', 'C', 'F', 'G', 'U', 'H', 'H', 'P', 'Q', 'O'],
    ['P', 'N', 'F', 'A', 'A', 'A', 'C', 'H', 'A', 'Z', 'K', 'J'],
    ['V', 'K', 'L', 'A', 'A', 'A', 'K', 'U', 'K', 'M', 'L', 'O'],
    ['N', 'G', 'R', 'A', 'A', 'A', 'E', 'E', 'H', 'F', 'Z', 'W'],
    ['I', 'A', 'D', 'Z', 'C', 'E', 'G', 'D', 'E', 'A', 'G', 'A'],
    ['H', 'H', 'H', 'T', 'N', 'L', 'S', 'V', 'F', 'W', 'T', 'V'],
    ['L', 'T', 'U', 'J', 'T', 'B', 'W', 'I', 'Q', 'Q', 'N', 'S']
]

result = find_squares(matrix)
print(result)
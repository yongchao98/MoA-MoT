def find_squares(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    min_side_length = 3
    result = []

    for i in range(rows):
        for j in range(cols):
            # Check for squares with side length >= 3
            for side_length in range(min_side_length, min(rows - i, cols - j) + 1):
                # Check if all characters in the square are the same
                char = matrix[i][j]
                is_square = True
                for x in range(i, i + side_length):
                    for y in range(j, j + side_length):
                        if matrix[x][y] != char:
                            is_square = False
                            break
                    if not is_square:
                        break
                if is_square:
                    result.append((i + side_length - 1, j + side_length - 1))
    
    # Since the answer is unique, we return the last found square's bottom-right corner
    return result[-1]

matrix = [
    ['K', 'G', 'H', 'D', 'K', 'T', 'M', 'T', 'X', 'T'],
    ['X', 'N', 'B', 'F', 'Q', 'L', 'R', 'J', 'Z', 'K'],
    ['D', 'T', 'Y', 'W', 'U', 'G', 'J', 'J', 'K', 'F'],
    ['C', 'L', 'X', 'M', 'W', 'X', 'U', 'J', 'E', 'B'],
    ['J', 'D', 'G', 'V', 'R', 'R', 'L', 'V', 'Z', 'N'],
    ['I', 'P', 'T', 'M', 'B', 'L', 'U', 'T', 'F', 'B'],
    ['Y', 'R', 'X', 'P', 'H', 'H', 'H', 'H', 'H', 'E'],
    ['N', 'X', 'M', 'N', 'H', 'H', 'H', 'H', 'H', 'C'],
    ['L', 'O', 'C', 'P', 'H', 'H', 'H', 'H', 'H', 'V'],
    ['N', 'C', 'L', 'X', 'H', 'H', 'H', 'H', 'H', 'L'],
    ['M', 'I', 'P', 'E', 'H', 'H', 'H', 'H', 'H', 'F'],
    ['K', 'L', 'C', 'Z', 'A', 'Y', 'G', 'V', 'W', 'Z']
]

result = find_squares(matrix)
print(result)
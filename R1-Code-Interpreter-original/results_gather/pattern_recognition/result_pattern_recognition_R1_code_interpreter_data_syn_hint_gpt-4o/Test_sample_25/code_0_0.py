def find_largest_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    max_size = 0
    bottom_right_corner = (-1, -1)

    for i in range(rows):
        for j in range(cols):
            # Check for potential squares starting from size 3x3
            for size in range(3, min(rows - i, cols - j) + 1):
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
                
                if is_square and size > max_size:
                    max_size = size
                    bottom_right_corner = (i + size - 1, j + size - 1)

    return bottom_right_corner

matrix = [
    ['N', 'H', 'M', 'K', 'H', 'L', 'X', 'G', 'S', 'R'],
    ['A', 'L', 'Z', 'D', 'Z', 'S', 'S', 'R', 'G', 'X'],
    ['N', 'O', 'P', 'I', 'I', 'I', 'I', 'B', 'I', 'H'],
    ['R', 'R', 'Y', 'I', 'I', 'I', 'I', 'H', 'D', 'V'],
    ['J', 'I', 'X', 'I', 'I', 'I', 'I', 'D', 'Y', 'I'],
    ['T', 'N', 'G', 'I', 'I', 'I', 'I', 'Z', 'L', 'J'],
    ['C', 'V', 'E', 'T', 'D', 'P', 'Q', 'T', 'Y', 'B'],
    ['C', 'J', 'A', 'N', 'S', 'S', 'J', 'S', 'O', 'E']
]

result = find_largest_square(matrix)
print(result)
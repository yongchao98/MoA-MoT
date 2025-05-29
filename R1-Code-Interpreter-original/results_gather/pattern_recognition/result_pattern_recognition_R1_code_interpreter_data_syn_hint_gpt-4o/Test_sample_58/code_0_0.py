def find_largest_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    max_size = 0
    bottom_right_corner = None

    for i in range(rows):
        for j in range(cols):
            # Check for squares starting from size 3
            for size in range(3, min(rows - i, cols - j) + 1):
                # Check if all characters in the square are the same
                char = matrix[i][j]
                valid_square = True
                for x in range(i, i + size):
                    for y in range(j, j + size):
                        if matrix[x][y] != char:
                            valid_square = False
                            break
                    if not valid_square:
                        break
                
                if valid_square and size > max_size:
                    max_size = size
                    bottom_right_corner = (i + size - 1, j + size - 1)

    return bottom_right_corner

matrix = [
    ['Z', 'A', 'W', 'G', 'J', 'I', 'L', 'W', 'H'],
    ['L', 'I', 'P', 'E', 'B', 'T', 'J', 'Z', 'C'],
    ['R', 'P', 'H', 'P', 'D', 'H', 'G', 'D', 'N'],
    ['J', 'J', 'J', 'D', 'N', 'T', 'P', 'N', 'G'],
    ['J', 'J', 'J', 'O', 'R', 'O', 'A', 'A', 'H'],
    ['J', 'J', 'J', 'M', 'C', 'V', 'P', 'T', 'J'],
    ['S', 'R', 'U', 'E', 'L', 'E', 'Z', 'K', 'F']
]

result = find_largest_square(matrix)
print(result)
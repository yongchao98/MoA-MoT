def find_largest_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    max_size = 0
    bottom_right_corner = None

    for i in range(rows):
        for j in range(cols):
            # Check for squares starting from size 3x3
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
    ['Q', 'N', 'M', 'K', 'W', 'Y', 'O', 'F'],
    ['B', 'U', 'C', 'Y', 'K', 'F', 'S', 'S'],
    ['I', 'M', 'P', 'M', 'K', 'Q', 'V', 'L'],
    ['Z', 'L', 'U', 'B', 'Q', 'D', 'I', 'Q'],
    ['E', 'Z', 'W', 'H', 'E', 'A', 'V', 'C'],
    ['A', 'P', 'Y', 'Y', 'Y', 'Y', 'A', 'R'],
    ['L', 'F', 'Y', 'Y', 'Y', 'Y', 'U', 'W'],
    ['P', 'F', 'Y', 'Y', 'Y', 'Y', 'G', 'C'],
    ['S', 'M', 'Y', 'Y', 'Y', 'Y', 'H', 'T'],
    ['Q', 'P', 'S', 'Z', 'C', 'K', 'K', 'Y'],
    ['L', 'L', 'O', 'N', 'O', 'G', 'A', 'A'],
    ['W', 'Z', 'G', 'S', 'F', 'Z', 'V', 'A'],
    ['P', 'X', 'P', 'P', 'X', 'N', 'N', 'X'],
    ['Z', 'G', 'Z', 'Q', 'K', 'S', 'Y', 'B']
]

result = find_largest_square(matrix)
print(result)
def find_largest_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    largest_square_bottom_right = None
    largest_size = 0

    for i in range(rows):
        for j in range(cols):
            # Check for squares starting from size 3x3
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
                
                if is_square and size > largest_size:
                    largest_size = size
                    largest_square_bottom_right = (i + size - 1, j + size - 1)

    return largest_square_bottom_right

matrix = [
    ['J', 'H', 'E', 'U', 'I', 'L', 'V', 'R', 'D'],
    ['G', 'A', 'S', 'N', 'G', 'N', 'P', 'M', 'K'],
    ['Z', 'H', 'Z', 'F', 'Z', 'L', 'R', 'B', 'C'],
    ['E', 'P', 'R', 'P', 'T', 'J', 'B', 'N', 'A'],
    ['S', 'U', 'A', 'L', 'X', 'U', 'A', 'P', 'V'],
    ['N', 'H', 'R', 'T', 'P', 'B', 'L', 'Y', 'M'],
    ['U', 'O', 'X', 'D', 'E', 'J', 'Q', 'X', 'E'],
    ['W', 'O', 'M', 'E', 'W', 'Q', 'H', 'J', 'H'],
    ['W', 'W', 'D', 'F', 'T', 'N', 'I', 'E', 'V'],
    ['C', 'I', 'Z', 'U', 'Y', 'L', 'M', 'Y', 'W'],
    ['Y', 'T', 'W', 'Q', 'K', 'T', 'G', 'O', 'N'],
    ['X', 'X', 'X', 'P', 'G', 'H', 'X', 'C', 'L'],
    ['X', 'X', 'X', 'K', 'A', 'X', 'A', 'K', 'T'],
    ['X', 'X', 'X', 'Y', 'I', 'H', 'I', 'B', 'T'],
    ['R', 'G', 'G', 'S', 'I', 'W', 'P', 'V', 'S']
]

result = find_largest_square(matrix)
print(result)
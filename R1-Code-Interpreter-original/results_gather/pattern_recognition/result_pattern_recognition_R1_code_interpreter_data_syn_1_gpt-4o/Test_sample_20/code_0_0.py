def find_largest_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    max_size = 0
    bottom_right_corner = (-1, -1)

    # Iterate over each cell as a potential top-left corner of a square
    for i in range(rows):
        for j in range(cols):
            # Try to form squares of increasing size
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
                
                # Update the largest square found
                if valid_square and size > max_size:
                    max_size = size
                    bottom_right_corner = (i + size - 1, j + size - 1)

    return bottom_right_corner

matrix = [
    ['N', 'B', 'M', 'H', 'I', 'U', 'W', 'F', 'Y', 'Q', 'S', 'D', 'Z'],
    ['L', 'I', 'Z', 'T', 'C', 'E', 'R', 'W', 'H', 'A', 'Y', 'E', 'Y'],
    ['H', 'P', 'C', 'R', 'E', 'B', 'T', 'Q', 'I', 'G', 'S', 'B', 'M'],
    ['N', 'F', 'E', 'W', 'T', 'Q', 'D', 'K', 'E', 'C', 'C', 'S', 'V'],
    ['G', 'M', 'A', 'S', 'Y', 'D', 'H', 'H', 'K', 'Z', 'E', 'Q', 'V'],
    ['I', 'B', 'D', 'A', 'V', 'S', 'B', 'M', 'V', 'O', 'W', 'F', 'M'],
    ['R', 'U', 'F', 'P', 'S', 'L', 'C', 'Q', 'F', 'Y', 'P', 'V', 'T'],
    ['G', 'C', 'T', 'U', 'T', 'P', 'N', 'A', 'J', 'P', 'P', 'E', 'T'],
    ['Q', 'J', 'K', 'R', 'T', 'J', 'V', 'D', 'K', 'B', 'I', 'K', 'L'],
    ['X', 'D', 'L', 'L', 'L', 'N', 'Z', 'Q', 'J', 'U', 'M', 'K', 'S'],
    ['M', 'R', 'L', 'L', 'L', 'M', 'N', 'Y', 'S', 'X', 'F', 'U', 'Z'],
    ['F', 'H', 'L', 'L', 'L', 'K', 'R', 'N', 'Q', 'E', 'Z', 'H', 'U'],
    ['C', 'D', 'J', 'E', 'Q', 'T', 'M', 'P', 'X', 'M', 'M', 'M', 'U']
]

result = find_largest_square(matrix)
print(result)
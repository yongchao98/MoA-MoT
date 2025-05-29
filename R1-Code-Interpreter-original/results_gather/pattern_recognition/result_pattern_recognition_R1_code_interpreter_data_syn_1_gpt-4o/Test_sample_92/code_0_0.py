def find_largest_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    largest_square_bottom_right = None
    largest_size = 0

    for i in range(rows):
        for j in range(cols):
            # Check for squares of increasing size starting from 3x3
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

result = find_largest_square(matrix)
print(result)
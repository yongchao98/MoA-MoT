def find_largest_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    max_size = 0
    bottom_right_corner = None

    for i in range(rows):
        for j in range(cols):
            # Check for squares of increasing size
            for size in range(3, min(rows - i, cols - j) + 1):
                # Check if all elements in the square are the same
                char = matrix[i][j]
                is_square = True
                for x in range(size):
                    for y in range(size):
                        if matrix[i + x][j + y] != char:
                            is_square = False
                            break
                    if not is_square:
                        break
                if is_square and size > max_size:
                    max_size = size
                    bottom_right_corner = (i + size - 1, j + size - 1)

    return bottom_right_corner

matrix = [
    ['J', 'P', 'N', 'C', 'G', 'Z', 'A', 'L', 'R', 'C', 'Y', 'J', 'F', 'O', 'V'],
    ['K', 'C', 'C', 'D', 'S', 'S', 'S', 'S', 'S', 'A', 'W', 'J', 'T', 'O', 'O'],
    ['H', 'M', 'G', 'Z', 'S', 'S', 'S', 'S', 'S', 'E', 'F', 'B', 'D', 'Q', 'Y'],
    ['M', 'V', 'F', 'X', 'S', 'S', 'S', 'S', 'S', 'U', 'W', 'U', 'Z', 'P', 'W'],
    ['B', 'M', 'W', 'R', 'S', 'S', 'S', 'S', 'S', 'L', 'K', 'V', 'U', 'M', 'R'],
    ['T', 'T', 'I', 'E', 'S', 'S', 'S', 'S', 'S', 'X', 'S', 'D', 'X', 'C', 'W'],
    ['L', 'X', 'B', 'T', 'Q', 'O', 'A', 'H', 'C', 'C', 'Y', 'F', 'G', 'P', 'Z'],
    ['R', 'B', 'G', 'Q', 'G', 'V', 'U', 'W', 'O', 'E', 'P', 'I', 'O', 'Y', 'C'],
    ['Z', 'K', 'C', 'H', 'U', 'M', 'M', 'P', 'Q', 'P', 'B', 'N', 'C', 'T', 'W'],
    ['H', 'J', 'F', 'Q', 'P', 'I', 'V', 'R', 'B', 'H', 'J', 'O', 'E', 'S', 'V'],
    ['R', 'H', 'T', 'I', 'J', 'N', 'J', 'J', 'N', 'N', 'G', 'X', 'S', 'H', 'G'],
    ['L', 'C', 'I', 'G', 'B', 'L', 'Q', 'Q', 'W', 'E', 'E', 'F', 'U', 'L', 'G'],
    ['I', 'Q', 'K', 'A', 'R', 'Z', 'J', 'B', 'X', 'K', 'U', 'A', 'R', 'Z', 'V'],
    ['Q', 'W', 'J', 'Q', 'J', 'A', 'D', 'N', 'A', 'D', 'W', 'I', 'F', 'U', 'O']
]

result = find_largest_square(matrix)
print(result)
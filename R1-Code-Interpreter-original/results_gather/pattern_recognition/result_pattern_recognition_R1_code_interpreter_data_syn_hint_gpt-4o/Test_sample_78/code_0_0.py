def find_largest_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    max_size = 0
    bottom_right_corner = (-1, -1)

    # Iterate over each cell as a potential top-left corner of a square
    for i in range(rows):
        for j in range(cols):
            # Check for squares of increasing size
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
                
                # Update the largest square found
                if is_square and size > max_size:
                    max_size = size
                    bottom_right_corner = (i + size - 1, j + size - 1)

    return bottom_right_corner

matrix = [
    ['V', 'J', 'Q', 'W', 'S', 'X', 'C', 'B'],
    ['T', 'X', 'R', 'K', 'K', 'Z', 'N', 'E'],
    ['J', 'T', 'X', 'R', 'A', 'S', 'L', 'X'],
    ['J', 'H', 'H', 'J', 'M', 'P', 'P', 'P'],
    ['Q', 'T', 'I', 'F', 'H', 'P', 'P', 'P'],
    ['A', 'K', 'T', 'X', 'R', 'P', 'P', 'P'],
    ['P', 'N', 'J', 'N', 'T', 'G', 'K', 'O']
]

result = find_largest_square(matrix)
print(result)
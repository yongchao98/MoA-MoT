def find_largest_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    max_size = 0
    bottom_right_corner = (-1, -1)

    # Iterate over each possible top-left corner of a square
    for i in range(rows):
        for j in range(cols):
            # Try to form squares of increasing size
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
                
                # If a valid square is found, update the largest square
                if is_square and size > max_size:
                    max_size = size
                    bottom_right_corner = (i + size - 1, j + size - 1)

    return bottom_right_corner

# Define the matrix
matrix = [
    ['G', 'S', 'J', 'W', 'E', 'M', 'Y', 'Y', 'C'],
    ['C', 'U', 'A', 'B', 'N', 'Z', 'Y', 'L', 'N'],
    ['E', 'L', 'F', 'F', 'F', 'F', 'Y', 'E', 'W'],
    ['Y', 'B', 'F', 'F', 'F', 'F', 'M', 'V', 'T'],
    ['H', 'P', 'F', 'F', 'F', 'F', 'R', 'N', 'H'],
    ['N', 'K', 'F', 'F', 'F', 'F', 'L', 'D', 'R'],
    ['Z', 'X', 'C', 'Y', 'I', 'W', 'D', 'V', 'W'],
    ['O', 'L', 'H', 'K', 'U', 'Y', 'K', 'U', 'P'],
    ['W', 'W', 'W', 'A', 'M', 'R', 'S', 'F', 'D'],
    ['C', 'C', 'R', 'C', 'A', 'P', 'L', 'A', 'X'],
    ['W', 'V', 'C', 'L', 'Q', 'T', 'O', 'J', 'X'],
    ['E', 'K', 'F', 'P', 'G', 'O', 'D', 'R', 'G'],
    ['G', 'A', 'I', 'B', 'I', 'A', 'M', 'A', 'C']
]

# Find and print the bottom-right corner of the largest square
print(find_largest_square(matrix))
def find_largest_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    max_size = 0
    bottom_right_corner = None

    # Iterate over each possible top-left corner of a square
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
                
                # If a valid square is found, update the largest size and bottom-right corner
                if valid_square and size > max_size:
                    max_size = size
                    bottom_right_corner = (i + size - 1, j + size - 1)

    return bottom_right_corner

# Define the matrix
matrix = [
    ['Z', 'J', 'B', 'B', 'O', 'W', 'C', 'O', 'S', 'V', 'L', 'K', 'A', 'M', 'K'],
    ['D', 'J', 'R', 'J', 'Q', 'T', 'G', 'J', 'U', 'E', 'N', 'M', 'R', 'A', 'U'],
    ['Y', 'Y', 'H', 'P', 'M', 'H', 'J', 'A', 'A', 'N', 'Y', 'W', 'P', 'A', 'M'],
    ['S', 'Z', 'Z', 'D', 'Y', 'W', 'J', 'H', 'J', 'S', 'M', 'X', 'B', 'V', 'R'],
    ['H', 'K', 'W', 'Q', 'V', 'U', 'A', 'L', 'B', 'A', 'I', 'D', 'E', 'V', 'P'],
    ['A', 'Y', 'M', 'C', 'N', 'X', 'X', 'X', 'U', 'S', 'R', 'P', 'X', 'O', 'O'],
    ['Z', 'Y', 'R', 'T', 'X', 'X', 'X', 'X', 'M', 'H', 'F', 'T', 'K', 'H', 'F'],
    ['P', 'I', 'Z', 'Z', 'Q', 'X', 'X', 'X', 'R', 'X', 'T', 'O', 'R', 'A', 'B'],
    ['M', 'E', 'V', 'R', 'X', 'B', 'U', 'K', 'W', 'L', 'Z', 'A', 'H', 'M', 'B'],
    ['P', 'S', 'I', 'X', 'F', 'E', 'Z', 'P', 'V', 'W', 'M', 'C', 'V', 'B', 'B'],
    ['E', 'K', 'U', 'A', 'P', 'B', 'M', 'T', 'P', 'J', 'T', 'A', 'S', 'L', 'A'],
    ['C', 'R', 'W', 'P', 'Q', 'R', 'Y', 'X', 'S', 'C', 'O', 'E', 'A', 'J', 'G'],
    ['T', 'S', 'H', 'U', 'U', 'I', 'G', 'J', 'L', 'U', 'B', 'Q', 'J', 'I', 'B']
]

# Find the bottom-right corner of the largest square
result = find_largest_square(matrix)
print(result)
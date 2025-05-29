def find_largest_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    max_square_size = 0
    bottom_right_corner = None

    # Iterate over each possible starting point for a square
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
                
                # If a valid square is found, update the largest square info
                if is_square and size > max_square_size:
                    max_square_size = size
                    bottom_right_corner = (i + size - 1, j + size - 1)

    return bottom_right_corner

# Define the matrix
matrix = [
    ['I', 'L', 'R', 'J', 'T', 'P', 'Z', 'C', 'D', 'J', 'B', 'L', 'H', 'X', 'M'],
    ['Q', 'N', 'D', 'D', 'G', 'T', 'L', 'R', 'Z', 'N', 'G', 'V', 'E', 'Y', 'T'],
    ['A', 'K', 'C', 'X', 'R', 'D', 'J', 'E', 'V', 'L', 'H', 'Y', 'A', 'U', 'P'],
    ['J', 'Z', 'R', 'L', 'V', 'U', 'H', 'G', 'Q', 'F', 'T', 'N', 'K', 'H', 'U'],
    ['I', 'Y', 'B', 'R', 'V', 'X', 'M', 'B', 'P', 'T', 'K', 'T', 'M', 'T', 'G'],
    ['N', 'W', 'X', 'Z', 'F', 'M', 'C', 'C', 'C', 'C', 'L', 'C', 'L', 'Y', 'V'],
    ['N', 'S', 'N', 'T', 'G', 'P', 'C', 'C', 'C', 'C', 'D', 'Q', 'D', 'G', 'F'],
    ['E', 'Z', 'V', 'Q', 'M', 'L', 'C', 'C', 'C', 'C', 'C', 'P', 'F', 'R', 'H'],
    ['C', 'K', 'C', 'O', 'P', 'U', 'C', 'C', 'C', 'C', 'K', 'D', 'H', 'Q', 'I'],
    ['R', 'H', 'B', 'H', 'B', 'H', 'W', 'M', 'O', 'F', 'T', 'J', 'C', 'A', 'E'],
    ['O', 'T', 'E', 'Q', 'W', 'Z', 'U', 'L', 'S', 'M', 'R', 'C', 'L', 'H', 'I'],
    ['O', 'E', 'C', 'S', 'Y', 'M', 'M', 'V', 'L', 'L', 'W', 'U', 'B', 'X', 'F'],
    ['L', 'Z', 'G', 'G', 'S', 'M', 'E', 'P', 'N', 'B', 'S', 'L', 'L', 'W', 'R'],
    ['Z', 'R', 'Z', 'G', 'R', 'M', 'Y', 'K', 'L', 'X', 'A', 'V', 'H', 'N', 'Z'],
    ['J', 'W', 'R', 'L', 'O', 'C', 'G', 'R', 'Z', 'X', 'S', 'F', 'B', 'Z', 'I']
]

# Find the bottom-right corner of the largest square
result = find_largest_square(matrix)
print(result)
def find_largest_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    max_side = 0
    bottom_right_corner = (0, 0)

    # Iterate over each possible starting point for a square
    for i in range(rows):
        for j in range(cols):
            # Try to form squares of increasing size
            for l in range(3, min(rows - i, cols - j) + 1):
                # Check if all characters in the square are the same
                char = matrix[i][j]
                valid_square = True
                for x in range(i, i + l):
                    for y in range(j, j + l):
                        if matrix[x][y] != char:
                            valid_square = False
                            break
                    if not valid_square:
                        break
                
                # If a valid square is found, check if it's the largest
                if valid_square and l > max_side:
                    max_side = l
                    bottom_right_corner = (i + l - 1, j + l - 1)

    return bottom_right_corner

# Define the matrix
matrix = [
    ['Y', 'V', 'F', 'S', 'F', 'A', 'T', 'P'],
    ['M', 'X', 'Z', 'P', 'Z', 'D', 'O', 'T'],
    ['U', 'O', 'O', 'X', 'U', 'M', 'E', 'Y'],
    ['E', 'T', 'A', 'G', 'Z', 'Q', 'A', 'L'],
    ['F', 'R', 'R', 'H', 'R', 'P', 'R', 'P'],
    ['G', 'T', 'K', 'T', 'Z', 'W', 'W', 'W'],
    ['M', 'S', 'E', 'L', 'Z', 'W', 'W', 'W'],
    ['E', 'G', 'Y', 'W', 'S', 'W', 'W', 'W'],
    ['Z', 'H', 'T', 'K', 'F', 'F', 'Y', 'S'],
    ['O', 'R', 'V', 'H', 'V', 'S', 'S', 'V'],
    ['F', 'J', 'X', 'J', 'X', 'X', 'B', 'R'],
    ['C', 'O', 'G', 'T', 'J', 'G', 'L', 'K'],
    ['A', 'A', 'R', 'R', 'D', 'L', 'R', 'I']
]

# Find the bottom-right corner of the largest square
result = find_largest_square(matrix)
print(result)
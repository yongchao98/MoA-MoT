def find_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    
    # We need at least a 3x3 matrix to form a square of side 3
    for i in range(rows - 2):
        for j in range(cols - 2):
            # Check for squares of increasing size starting from 3x3
            for size in range(3, min(rows - i, cols - j) + 1):
                char = matrix[i][j]
                is_square = True
                # Check if all elements in the square are the same
                for x in range(i, i + size):
                    for y in range(j, j + size):
                        if matrix[x][y] != char:
                            is_square = False
                            break
                    if not is_square:
                        break
                if is_square:
                    # Record the bottom-right corner of the square
                    bottom_right = (i + size - 1, j + size - 1)
                    return bottom_right

# Define the matrix
matrix = [
    ['G', 'Y', 'A', 'D', 'H', 'G', 'C', 'P', 'H', 'F', 'E', 'E', 'F', 'G', 'F'],
    ['X', 'K', 'G', 'F', 'T', 'E', 'D', 'D', 'E', 'Y', 'H', 'N', 'J', 'R', 'P'],
    ['V', 'L', 'B', 'F', 'A', 'I', 'I', 'I', 'I', 'H', 'O', 'N', 'O', 'R', 'H'],
    ['E', 'N', 'G', 'U', 'Z', 'I', 'I', 'I', 'I', 'W', 'J', 'O', 'M', 'V', 'F'],
    ['J', 'E', 'M', 'F', 'U', 'I', 'I', 'I', 'I', 'A', 'Y', 'B', 'T', 'F', 'E'],
    ['D', 'R', 'A', 'H', 'R', 'I', 'I', 'I', 'I', 'X', 'R', 'M', 'I', 'S', 'R'],
    ['O', 'O', 'O', 'C', 'X', 'X', 'O', 'R', 'F', 'F', 'K', 'L', 'Y', 'G', 'N'],
    ['P', 'R', 'X', 'Y', 'M', 'G', 'W', 'L', 'G', 'A', 'R', 'B', 'G', 'I', 'F']
]

# Find the bottom-right corner of the square
result = find_square(matrix)
print(result)
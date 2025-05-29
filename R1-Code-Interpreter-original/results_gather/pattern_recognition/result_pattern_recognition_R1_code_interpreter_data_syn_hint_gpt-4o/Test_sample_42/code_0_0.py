def find_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    
    # Minimum side length of the square
    min_side = 3
    
    # Iterate over each possible starting point for a square
    for i in range(rows):
        for j in range(cols):
            # Check possible square sizes starting from (i, j)
            for side in range(min_side, min(rows - i, cols - j) + 1):
                # Check if all characters in the square are the same
                char = matrix[i][j]
                is_square = True
                for x in range(i, i + side):
                    for y in range(j, j + side):
                        if matrix[x][y] != char:
                            is_square = False
                            break
                    if not is_square:
                        break
                
                # If a valid square is found, record the bottom-right corner
                if is_square:
                    bottom_right = (i + side - 1, j + side - 1)
                    return bottom_right

# Define the matrix
matrix = [
    ['W', 'E', 'B', 'P', 'M', 'L'],
    ['E', 'R', 'Q', 'T', 'S', 'L'],
    ['Z', 'V', 'V', 'V', 'D', 'H'],
    ['R', 'V', 'V', 'V', 'F', 'A'],
    ['K', 'V', 'V', 'V', 'N', 'F'],
    ['R', 'K', 'R', 'V', 'D', 'J'],
    ['J', 'W', 'T', 'E', 'W', 'J'],
    ['E', 'A', 'M', 'R', 'E', 'K']
]

# Find the bottom-right corner of the square
result = find_square(matrix)
print(result)
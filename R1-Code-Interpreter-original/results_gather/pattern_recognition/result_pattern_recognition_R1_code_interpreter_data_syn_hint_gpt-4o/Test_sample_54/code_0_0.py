def find_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    
    # Minimum side length of the square
    min_side = 3
    
    # Iterate over each possible top-left corner of the square
    for i in range(rows):
        for j in range(cols):
            # Check possible square sizes starting from 3x3
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
                
                # If a valid square is found, return the bottom-right corner
                if is_square:
                    bottom_right = (i + side - 1, j + side - 1)
                    return bottom_right

# Define the matrix
matrix = [
    ['J', 'R', 'Z', 'T', 'F', 'B', 'R'],
    ['N', 'F', 'J', 'S', 'Q', 'U', 'D'],
    ['F', 'B', 'H', 'U', 'X', 'L', 'W'],
    ['Y', 'M', 'K', 'U', 'L', 'H', 'N'],
    ['C', 'U', 'H', 'I', 'J', 'B', 'Q'],
    ['I', 'N', 'N', 'N', 'T', 'Z', 'I'],
    ['A', 'N', 'N', 'N', 'E', 'D', 'H'],
    ['B', 'N', 'N', 'N', 'X', 'U', 'R'],
    ['O', 'X', 'O', 'P', 'Q', 'R', 'I']
]

# Find the bottom-right corner of the square
result = find_square(matrix)
print(result)
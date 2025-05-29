def find_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    
    # Minimum side length of the square
    min_side = 3
    
    # Iterate over each possible starting point for a square
    for r in range(rows):
        for c in range(cols):
            # Check possible square sizes
            for side in range(min_side, min(rows - r, cols - c) + 1):
                # Check if all characters in the square are the same
                char = matrix[r][c]
                is_square = True
                for i in range(r, r + side):
                    for j in range(c, c + side):
                        if matrix[i][j] != char:
                            is_square = False
                            break
                    if not is_square:
                        break
                
                # If a valid square is found, return the bottom-right corner
                if is_square:
                    return [r + side - 1, c + side - 1]

# Define the matrix
matrix = [
    ['S', 'L', 'R', 'S', 'D', 'A', 'J', 'D', 'K', 'K'],
    ['D', 'Z', 'O', 'T', 'M', 'N', 'O', 'I', 'I', 'L'],
    ['T', 'J', 'X', 'K', 'F', 'F', 'F', 'F', 'F', 'W'],
    ['Q', 'Z', 'A', 'A', 'F', 'F', 'F', 'F', 'F', 'Q'],
    ['L', 'W', 'R', 'U', 'F', 'F', 'F', 'F', 'F', 'S'],
    ['Y', 'M', 'X', 'M', 'F', 'F', 'F', 'F', 'F', 'N'],
    ['P', 'E', 'I', 'E', 'F', 'F', 'F', 'F', 'F', 'F'],
    ['X', 'F', 'K', 'C', 'W', 'J', 'R', 'W', 'N', 'P'],
    ['S', 'B', 'P', 'M', 'V', 'J', 'S', 'S', 'L', 'L'],
    ['C', 'S', 'M', 'U', 'Z', 'U', 'U', 'T', 'F', 'W'],
    ['E', 'R', 'J', 'A', 'S', 'P', 'T', 'Z', 'L', 'S'],
    ['X', 'C', 'L', 'R', 'Y', 'L', 'B', 'X', 'W', 'I']
]

# Find and print the bottom-right corner of the square
print(find_square(matrix))
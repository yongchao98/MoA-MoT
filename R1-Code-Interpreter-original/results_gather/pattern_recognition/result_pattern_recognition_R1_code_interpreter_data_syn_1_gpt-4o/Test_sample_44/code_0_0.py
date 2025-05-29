def find_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    
    # Iterate over each possible top-left corner of a square
    for i in range(rows):
        for j in range(cols):
            # Check for squares of size 3x3 or larger
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
                
                # If a valid square is found, record the bottom-right corner
                if valid_square:
                    bottom_right = (i + size - 1, j + size - 1)
                    return bottom_right

# Define the matrix
matrix = [
    ['S', 'M', 'Q', 'A', 'J', 'N', 'Z', 'I', 'B', 'F', 'H', 'W', 'M', 'W'],
    ['A', 'V', 'O', 'A', 'Z', 'Z', 'Z', 'M', 'X', 'B', 'C', 'N', 'M', 'V'],
    ['X', 'Z', 'Z', 'D', 'Z', 'Z', 'Z', 'S', 'V', 'J', 'N', 'R', 'D', 'L'],
    ['S', 'Z', 'T', 'O', 'Z', 'Z', 'Z', 'P', 'M', 'K', 'K', 'K', 'D', 'T'],
    ['B', 'W', 'O', 'Q', 'E', 'I', 'B', 'A', 'K', 'A', 'Q', 'X', 'N', 'T'],
    ['K', 'O', 'V', 'B', 'U', 'G', 'O', 'B', 'Q', 'S', 'F', 'R', 'U', 'J']
]

# Find and print the bottom-right corner of the square
print(find_square(matrix))
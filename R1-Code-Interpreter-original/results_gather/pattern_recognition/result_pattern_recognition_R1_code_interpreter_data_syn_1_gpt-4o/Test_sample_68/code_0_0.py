def find_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    
    # Iterate over each possible top-left corner of a square
    for i in range(rows):
        for j in range(cols):
            # Check for squares of increasing size starting from 3x3
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
                
                # If a valid square is found, return the bottom-right corner
                if valid_square:
                    return [i + size - 1, j + size - 1]

# Define the matrix
matrix = [
    ['T', 'X', 'A', 'L', 'J', 'V', 'X', 'Y', 'Y', 'L', 'R', 'U', 'X', 'I'],
    ['E', 'G', 'D', 'U', 'Y', 'T', 'C', 'F', 'F', 'F', 'F', 'L', 'I', 'U'],
    ['K', 'X', 'S', 'B', 'E', 'R', 'D', 'F', 'F', 'F', 'F', 'U', 'S', 'R'],
    ['F', 'G', 'U', 'K', 'P', 'F', 'O', 'F', 'F', 'F', 'F', 'F', 'C', 'E'],
    ['C', 'K', 'W', 'Z', 'C', 'V', 'Y', 'F', 'F', 'F', 'F', 'F', 'K', 'D'],
    ['H', 'E', 'N', 'J', 'Y', 'T', 'O', 'V', 'M', 'U', 'O', 'Y', 'J', 'J'],
    ['R', 'A', 'Y', 'Y', 'N', 'K', 'K', 'I', 'K', 'I', 'O', 'P', 'D', 'P'],
    ['E', 'M', 'G', 'B', 'E', 'E', 'G', 'X', 'G', 'P', 'T', 'R', 'K', 'K']
]

# Find and print the bottom-right corner of the square
print(find_square(matrix))
def find_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    
    # Iterate over each possible top-left corner of a square
    for i in range(rows - 2):
        for j in range(cols - 2):
            # Check for squares of increasing size starting from 3x3
            for size in range(3, min(rows - i, cols - j) + 1):
                # Check if all characters in the current square are the same
                char = matrix[i][j]
                is_square = True
                for x in range(i, i + size):
                    for y in range(j, j + size):
                        if matrix[x][y] != char:
                            is_square = False
                            break
                    if not is_square:
                        break
                
                # If a valid square is found, return the bottom-right corner
                if is_square:
                    return [i + size - 1, j + size - 1]

# Define the matrix
matrix = [
    ['F', 'P', 'F', 'Y', 'F', 'H', 'X', 'U', 'E', 'W', 'R', 'Z', 'C'],
    ['O', 'M', 'G', 'M', 'J', 'E', 'A', 'G', 'R', 'B', 'C', 'Z', 'U'],
    ['K', 'Z', 'Q', 'Z', 'K', 'F', 'F', 'X', 'Q', 'S', 'B', 'B', 'B'],
    ['D', 'W', 'W', 'A', 'G', 'U', 'U', 'M', 'X', 'Q', 'B', 'B', 'B'],
    ['Y', 'W', 'C', 'F', 'Q', 'F', 'X', 'Q', 'W', 'X', 'B', 'B', 'B'],
    ['R', 'N', 'V', 'R', 'P', 'K', 'Q', 'Q', 'S', 'W', 'M', 'W', 'I'],
    ['G', 'Q', 'N', 'A', 'E', 'H', 'R', 'B', 'O', 'U', 'N', 'N', 'J']
]

# Find and print the bottom-right corner of the square
print(find_square(matrix))
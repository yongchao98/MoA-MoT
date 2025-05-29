def find_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    
    # Iterate over each possible top-left corner of a square
    for i in range(rows - 2):
        for j in range(cols - 2):
            # Check for squares of increasing size
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
                if is_square:
                    # Return the bottom-right corner of the square
                    return [i + size - 1, j + size - 1]

# Define the matrix
matrix = [
    ['P', 'Y', 'L', 'I', 'G', 'H', 'C', 'G', 'I', 'W', 'C', 'C', 'Q'],
    ['V', 'Y', 'Q', 'Q', 'Q', 'Q', 'Q', 'K', 'Q', 'Q', 'O', 'Z', 'H'],
    ['K', 'S', 'Q', 'Q', 'Q', 'Q', 'Q', 'F', 'D', 'D', 'J', 'N', 'T'],
    ['B', 'P', 'Q', 'Q', 'Q', 'Q', 'Q', 'W', 'R', 'I', 'B', 'T', 'I'],
    ['Y', 'L', 'Q', 'Q', 'Q', 'Q', 'Q', 'Y', 'U', 'Q', 'S', 'Z', 'V'],
    ['D', 'F', 'Q', 'Q', 'Q', 'Q', 'Q', 'B', 'O', 'E', 'R', 'R', 'U'],
    ['Y', 'N', 'F', 'U', 'P', 'B', 'J', 'M', 'Y', 'E', 'U', 'H', 'X'],
    ['A', 'H', 'Z', 'O', 'C', 'B', 'I', 'I', 'D', 'D', 'W', 'X', 'F'],
    ['D', 'E', 'A', 'F', 'G', 'W', 'W', 'K', 'J', 'E', 'A', 'T', 'U'],
    ['R', 'Q', 'F', 'B', 'H', 'F', 'M', 'M', 'K', 'N', 'F', 'O', 'N']
]

# Find and print the bottom-right corner of the square
print(find_square(matrix))
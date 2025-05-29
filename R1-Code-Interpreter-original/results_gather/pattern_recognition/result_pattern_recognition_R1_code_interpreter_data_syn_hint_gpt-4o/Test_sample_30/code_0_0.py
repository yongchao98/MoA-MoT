def find_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    
    for i in range(rows - 2):
        for j in range(cols - 2):
            char = matrix[i][j]
            # Check for squares of increasing size starting from 3
            for size in range(3, min(rows - i, cols - j) + 1):
                is_square = True
                # Check the square of current size
                for x in range(size):
                    for y in range(size):
                        if matrix[i + x][j + y] != char:
                            is_square = False
                            break
                    if not is_square:
                        break
                if is_square:
                    # Return the bottom-right corner of the square
                    return [i + size - 1, j + size - 1]

# Define the matrix
matrix = [
    ['Q', 'Q', 'N', 'P', 'I', 'O'],
    ['N', 'M', 'P', 'W', 'M', 'K'],
    ['F', 'A', 'W', 'M', 'W', 'S'],
    ['C', 'E', 'E', 'E', 'D', 'V'],
    ['J', 'E', 'E', 'E', 'F', 'Z'],
    ['W', 'E', 'E', 'E', 'N', 'I'],
    ['I', 'Y', 'K', 'L', 'X', 'J'],
    ['P', 'W', 'L', 'G', 'H', 'P'],
    ['Q', 'U', 'M', 'G', 'D', 'Q'],
    ['U', 'M', 'W', 'H', 'D', 'U'],
    ['F', 'U', 'M', 'X', 'A', 'S']
]

# Find and print the bottom-right corner of the square
print(find_square(matrix))
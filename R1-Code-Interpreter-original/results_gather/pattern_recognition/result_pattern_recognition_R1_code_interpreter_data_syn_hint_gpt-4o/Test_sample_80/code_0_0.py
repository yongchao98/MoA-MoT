def find_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    
    # Iterate over each possible starting point for a 3x3 square
    for i in range(rows - 2):
        for j in range(cols - 2):
            # Check for squares of increasing size starting from 3x3
            max_size = min(rows - i, cols - j)
            for size in range(3, max_size + 1):
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
                    # Record the bottom-right corner of the square
                    bottom_right = (i + size - 1, j + size - 1)
                    return bottom_right

# Define the matrix
matrix = [
    ['Q', 'Y', 'Q', 'R', 'K', 'V', 'A', 'Q', 'Z', 'C', 'J'],
    ['S', 'K', 'F', 'C', 'N', 'U', 'X', 'O', 'X', 'X', 'I'],
    ['M', 'O', 'B', 'Q', 'Q', 'P', 'G', 'B', 'I', 'Y', 'W'],
    ['I', 'Y', 'J', 'B', 'G', 'P', 'B', 'D', 'Z', 'G', 'G'],
    ['I', 'Z', 'U', 'L', 'Y', 'P', 'M', 'M', 'M', 'M', 'T'],
    ['C', 'V', 'P', 'H', 'M', 'V', 'M', 'M', 'M', 'M', 'P'],
    ['M', 'H', 'I', 'T', 'Q', 'U', 'M', 'M', 'M', 'M', 'U'],
    ['Z', 'D', 'L', 'X', 'G', 'D', 'M', 'M', 'M', 'M', 'D'],
    ['H', 'Q', 'W', 'W', 'P', 'N', 'U', 'X', 'I', 'L', 'W'],
    ['L', 'G', 'I', 'Q', 'Y', 'P', 'E', 'W', 'B', 'U', 'Y'],
    ['B', 'Z', 'F', 'S', 'K', 'E', 'V', 'K', 'D', 'T', 'M'],
    ['B', 'Z', 'P', 'V', 'J', 'U', 'U', 'Y', 'H', 'V', 'Q']
]

# Find and print the bottom-right corner of the square
print(find_square(matrix))
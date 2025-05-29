def find_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    
    # Function to check if a square of size `size` can be formed starting at (r, c)
    def is_square(r, c, size):
        char = matrix[r][c]
        for i in range(size):
            for j in range(size):
                if matrix[r + i][c + j] != char:
                    return False
        return True
    
    # Iterate over each possible top-left corner of a square
    for r in range(rows):
        for c in range(cols):
            # Check for squares of increasing size
            for size in range(3, min(rows - r, cols - c) + 1):
                if is_square(r, c, size):
                    # Return the bottom-right corner of the first found square
                    return [r + size - 1, c + size - 1]

# Define the matrix
matrix = [
    ['R', 'L', 'L', 'L', 'Q', 'A', 'K', 'G', 'G', 'V'],
    ['H', 'L', 'L', 'L', 'R', 'U', 'Z', 'X', 'P', 'H'],
    ['N', 'L', 'L', 'L', 'N', 'J', 'F', 'V', 'P', 'J'],
    ['E', 'V', 'P', 'B', 'O', 'P', 'H', 'H', 'T', 'A'],
    ['G', 'R', 'Q', 'Y', 'K', 'V', 'Y', 'Y', 'H', 'Q'],
    ['Y', 'G', 'R', 'N', 'K', 'Y', 'R', 'H', 'Y', 'E'],
    ['J', 'E', 'M', 'A', 'S', 'L', 'F', 'I', 'Q', 'N']
]

# Find and print the bottom-right corner of the square
print(find_square(matrix))
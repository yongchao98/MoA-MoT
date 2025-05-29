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
            size = 3
            while r + size <= rows and c + size <= cols:
                if is_square(r, c, size):
                    # If a square is found, record the bottom-right corner
                    bottom_right = (r + size - 1, c + size - 1)
                size += 1
    
    return bottom_right

matrix = [
    ['K', 'X', 'I', 'R', 'T', 'M', 'J', 'C', 'W', 'A', 'Z'],
    ['G', 'U', 'P', 'M', 'T', 'U', 'E', 'V', 'G', 'I', 'M'],
    ['E', 'C', 'L', 'V', 'J', 'J', 'J', 'F', 'T', 'K', 'L'],
    ['A', 'J', 'M', 'W', 'J', 'J', 'J', 'Q', 'T', 'O', 'P'],
    ['C', 'E', 'U', 'X', 'J', 'J', 'J', 'X', 'A', 'U', 'G'],
    ['O', 'L', 'H', 'F', 'F', 'O', 'Y', 'Y', 'Z', 'O', 'K'],
    ['C', 'K', 'C', 'T', 'Z', 'G', 'E', 'G', 'M', 'L', 'D'],
    ['P', 'G', 'H', 'Y', 'U', 'C', 'U', 'Y', 'S', 'E', 'Y'],
    ['R', 'U', 'Z', 'H', 'R', 'C', 'E', 'Q', 'R', 'J', 'V'],
    ['R', 'Q', 'Y', 'T', 'K', 'X', 'I', 'V', 'U', 'I', 'S'],
    ['D', 'A', 'R', 'H', 'Q', 'N', 'H', 'A', 'B', 'I', 'Y'],
    ['T', 'K', 'X', 'F', 'D', 'D', 'G', 'A', 'T', 'M', 'I']
]

result = find_square(matrix)
print(result)
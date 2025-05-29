def find_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    
    # Start checking for squares of size 3x3 and larger
    for size in range(3, min(rows, cols) + 1):
        for i in range(rows - size + 1):
            for j in range(cols - size + 1):
                # Check if all elements in the square are the same
                char = matrix[i][j]
                is_square = True
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
    ['C', 'L', 'F', 'T', 'H', 'P', 'R', 'T', 'F', 'U', 'W', 'B'],
    ['R', 'A', 'Q', 'O', 'U', 'H', 'Z', 'I', 'H', 'S', 'E', 'D'],
    ['K', 'A', 'Z', 'O', 'Z', 'Q', 'K', 'J', 'I', 'I', 'I', 'D'],
    ['L', 'O', 'Q', 'X', 'M', 'J', 'S', 'I', 'I', 'I', 'I', 'G'],
    ['R', 'C', 'U', 'J', 'A', 'I', 'D', 'Q', 'I', 'I', 'I', 'H'],
    ['G', 'D', 'J', 'E', 'P', 'L', 'H', 'L', 'Q', 'O', 'C', 'U']
]

# Find and print the bottom-right corner of the square
print(find_square(matrix))
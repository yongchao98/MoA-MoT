def find_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    
    # We need to find the largest square with the same character
    for size in range(min(rows, cols), 2, -1):  # Start from the largest possible size
        for i in range(rows - size + 1):
            for j in range(cols - size + 1):
                # Check if all characters in the square are the same
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
    ['U', 'L', 'M', 'H', 'P', 'C', 'C', 'S', 'T', 'Y', 'L', 'P'],
    ['A', 'C', 'L', 'O', 'U', 'R', 'G', 'G', 'X', 'I', 'O', 'T'],
    ['H', 'F', 'W', 'S', 'C', 'B', 'C', 'W', 'D', 'V', 'U', 'G'],
    ['K', 'M', 'M', 'M', 'M', 'Q', 'S', 'Y', 'Z', 'U', 'Q', 'E'],
    ['C', 'M', 'M', 'M', 'M', 'O', 'O', 'E', 'S', 'A', 'P', 'I'],
    ['N', 'M', 'M', 'M', 'M', 'F', 'W', 'V', 'Z', 'M', 'D', 'Q'],
    ['Q', 'M', 'M', 'M', 'M', 'S', 'W', 'R', 'Y', 'T', 'O', 'N'],
    ['N', 'Z', 'C', 'Z', 'Z', 'V', 'N', 'C', 'R', 'L', 'F', 'P']
]

# Find and print the bottom-right corner of the largest square
print(find_square(matrix))
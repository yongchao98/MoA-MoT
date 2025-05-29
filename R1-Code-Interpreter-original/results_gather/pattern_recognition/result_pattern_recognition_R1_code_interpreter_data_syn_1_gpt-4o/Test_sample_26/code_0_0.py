def find_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    
    # We need at least a 3x3 matrix to form a square of side 3
    for size in range(3, min(rows, cols) + 1):
        for i in range(rows - size + 1):
            for j in range(cols - size + 1):
                # Check if all elements in the square are the same
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
    ['L', 'H', 'P', 'T', 'T', 'K', 'A', 'S'],
    ['B', 'W', 'W', 'W', 'D', 'I', 'J', 'W'],
    ['U', 'W', 'W', 'W', 'A', 'A', 'E', 'A'],
    ['O', 'W', 'W', 'W', 'T', 'T', 'Y', 'S'],
    ['D', 'M', 'Z', 'V', 'D', 'Z', 'Y', 'V'],
    ['S', 'C', 'T', 'W', 'H', 'C', 'Z', 'P'],
    ['K', 'M', 'N', 'Q', 'E', 'H', 'D', 'N'],
    ['E', 'K', 'S', 'X', 'A', 'U', 'I', 'D'],
    ['F', 'F', 'D', 'S', 'V', 'J', 'P', 'J'],
    ['I', 'I', 'S', 'I', 'T', 'I', 'V', 'Z'],
    ['R', 'Z', 'L', 'W', 'Y', 'J', 'D', 'E']
]

# Find the bottom-right corner of the square
result = find_square(matrix)
print(result)
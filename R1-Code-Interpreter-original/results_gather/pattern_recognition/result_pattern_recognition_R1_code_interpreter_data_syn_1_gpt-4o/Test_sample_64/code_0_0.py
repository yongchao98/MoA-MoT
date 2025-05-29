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
    ['F', 'N', 'U', 'V', 'P', 'Q'],
    ['C', 'I', 'I', 'I', 'V', 'D'],
    ['Y', 'I', 'I', 'I', 'M', 'R'],
    ['Z', 'I', 'I', 'I', 'S', 'J'],
    ['C', 'N', 'C', 'H', 'L', 'O'],
    ['O', 'S', 'A', 'Z', 'Y', 'H'],
    ['C', 'E', 'J', 'R', 'Y', 'Z'],
    ['X', 'A', 'V', 'N', 'P', 'E'],
    ['O', 'K', 'D', 'Q', 'Q', 'I'],
    ['E', 'H', 'G', 'N', 'P', 'F']
]

# Find and print the bottom-right corner of the square
print(find_square(matrix))
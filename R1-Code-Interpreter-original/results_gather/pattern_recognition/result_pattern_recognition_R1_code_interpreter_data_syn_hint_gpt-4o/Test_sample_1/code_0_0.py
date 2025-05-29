def find_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    
    for i in range(rows):
        for j in range(cols):
            # Check for squares of size 3x3 or larger
            for size in range(3, min(rows - i, cols - j) + 1):
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
    ['M', 'F', 'Q', 'G', 'Q', 'B', 'D', 'Z', 'X', 'D', 'N', 'J'],
    ['W', 'W', 'I', 'D', 'M', 'O', 'B', 'U', 'S', 'S', 'S', 'Q'],
    ['F', 'N', 'V', 'I', 'F', 'S', 'T', 'P', 'S', 'S', 'S', 'V'],
    ['V', 'E', 'Y', 'B', 'Y', 'V', 'K', 'H', 'S', 'S', 'S', 'H'],
    ['I', 'Q', 'D', 'W', 'Q', 'L', 'D', 'R', 'P', 'B', 'D', 'O'],
    ['Y', 'C', 'G', 'D', 'E', 'G', 'O', 'N', 'I', 'D', 'D', 'R'],
    ['R', 'Q', 'K', 'R', 'C', 'I', 'W', 'G', 'S', 'P', 'A', 'H']
]

# Find and print the bottom-right corner of the square
print(find_square(matrix))
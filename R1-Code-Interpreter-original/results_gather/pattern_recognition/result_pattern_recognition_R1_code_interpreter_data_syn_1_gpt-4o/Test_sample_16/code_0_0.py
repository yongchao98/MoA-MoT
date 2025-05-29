def find_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    
    # Iterate over each possible top-left corner of a square
    for i in range(rows - 2):
        for j in range(cols - 2):
            # Check for squares of size 3x3 or larger
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
                
                # If a valid square is found, return the bottom-right corner
                if is_square:
                    return [i + size - 1, j + size - 1]

# Define the matrix
matrix = [
    ['E', 'O', 'T', 'I', 'A', 'S'],
    ['D', 'Y', 'R', 'B', 'U', 'O'],
    ['Y', 'M', 'S', 'Z', 'Z', 'Z'],
    ['N', 'G', 'H', 'Z', 'Z', 'Z'],
    ['Q', 'S', 'X', 'Z', 'Z', 'Z'],
    ['V', 'P', 'M', 'K', 'Z', 'J'],
    ['E', 'X', 'X', 'M', 'F', 'N']
]

# Find and print the bottom-right corner of the square
print(find_square(matrix))
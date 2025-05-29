def find_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    
    # We need at least a 3x3 area to form a square of side 3
    for i in range(rows - 2):
        for j in range(cols - 2):
            # Check for squares of increasing size starting from 3
            for size in range(3, min(rows - i, cols - j) + 1):
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
    ['I', 'E', 'L', 'F', 'L', 'B'],
    ['Q', 'V', 'L', 'Q', 'F', 'X'],
    ['T', 'P', 'K', 'J', 'U', 'R'],
    ['D', 'Q', 'D', 'W', 'A', 'Y'],
    ['W', 'S', 'S', 'S', 'S', 'S'],
    ['A', 'S', 'S', 'S', 'Y', 'E'],
    ['B', 'S', 'S', 'S', 'G', 'J'],
    ['V', 'R', 'T', 'I', 'X', 'M']
]

# Find and print the bottom-right corner of the square
print(find_square(matrix))
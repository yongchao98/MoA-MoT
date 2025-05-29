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
    ['O', 'P', 'R', 'F', 'A', 'L', 'M', 'G'],
    ['C', 'Y', 'D', 'Z', 'F', 'V', 'G', 'V'],
    ['Q', 'I', 'J', 'B', 'I', 'V', 'B', 'E'],
    ['F', 'A', 'F', 'V', 'V', 'V', 'V', 'L'],
    ['D', 'B', 'Z', 'V', 'V', 'V', 'V', 'P'],
    ['Y', 'I', 'O', 'V', 'V', 'V', 'V', 'T'],
    ['K', 'K', 'J', 'V', 'V', 'V', 'V', 'Q'],
    ['C', 'V', 'O', 'F', 'E', 'Z', 'J', 'C'],
    ['J', 'F', 'D', 'V', 'C', 'F', 'S', 'A']
]

# Find and print the bottom-right corner of the square
print(find_square(matrix))
def find_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    
    # Iterate over each element in the matrix
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
    ['M', 'U', 'Z', 'H', 'H', 'H', 'H', 'N', 'T'],
    ['A', 'H', 'O', 'H', 'H', 'H', 'H', 'M', 'N'],
    ['R', 'A', 'J', 'H', 'H', 'H', 'H', 'Q', 'G'],
    ['K', 'J', 'X', 'H', 'H', 'H', 'H', 'V', 'L'],
    ['V', 'Y', 'L', 'Q', 'E', 'K', 'D', 'P', 'Y'],
    ['O', 'J', 'P', 'E', 'X', 'I', 'P', 'V', 'R'],
    ['T', 'H', 'V', 'B', 'Z', 'V', 'P', 'A', 'A'],
    ['X', 'D', 'C', 'H', 'S', 'D', 'H', 'S', 'C'],
    ['P', 'R', 'Q', 'Y', 'Q', 'O', 'Y', 'N', 'C'],
    ['X', 'A', 'O', 'M', 'R', 'K', 'D', 'X', 'I'],
    ['C', 'A', 'H', 'X', 'G', 'L', 'S', 'L', 'F'],
    ['Q', 'O', 'R', 'U', 'P', 'N', 'B', 'I', 'X'],
    ['D', 'R', 'Y', 'S', 'W', 'Q', 'W', 'W', 'E'],
    ['Q', 'V', 'Y', 'Y', 'A', 'Q', 'V', 'H', 'O'],
    ['N', 'I', 'V', 'N', 'V', 'X', 'O', 'Q', 'D']
]

# Find and print the bottom-right corner of the square
print(find_square(matrix))
def find_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    
    for i in range(rows):
        for j in range(cols):
            # Check for squares of increasing size starting from 3
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
    ['U', 'D', 'L', 'G', 'V', 'L', 'V'],
    ['A', 'K', 'L', 'K', 'E', 'N', 'N'],
    ['L', 'Z', 'O', 'A', 'V', 'A', 'M'],
    ['M', 'M', 'M', 'Z', 'K', 'W', 'S'],
    ['M', 'M', 'M', 'J', 'C', 'M', 'Y'],
    ['M', 'M', 'M', 'C', 'W', 'Q', 'D'],
    ['M', 'Y', 'Q', 'P', 'I', 'X', 'Z'],
    ['A', 'F', 'F', 'U', 'D', 'U', 'G'],
    ['R', 'L', 'N', 'P', 'P', 'U', 'K'],
    ['P', 'C', 'M', 'C', 'R', 'N', 'E'],
    ['V', 'E', 'O', 'J', 'B', 'R', 'I'],
    ['X', 'Y', 'C', 'R', 'N', 'R', 'R'],
    ['U', 'C', 'F', 'Z', 'Q', 'T', 'A'],
    ['Z', 'W', 'H', 'X', 'P', 'M', 'P']
]

# Find and print the bottom-right corner of the square
print(find_square(matrix))
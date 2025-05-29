def find_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    
    # Iterate over each possible top-left corner of a square
    for i in range(rows - 2):
        for j in range(cols - 2):
            # Check for squares of increasing size
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
    ['J', 'P', 'N', 'C', 'G', 'Z', 'A', 'L', 'R', 'C', 'Y', 'J', 'F', 'O', 'V'],
    ['K', 'C', 'C', 'D', 'S', 'S', 'S', 'S', 'S', 'A', 'W', 'J', 'T', 'O', 'O'],
    ['H', 'M', 'G', 'Z', 'S', 'S', 'S', 'S', 'S', 'E', 'F', 'B', 'D', 'Q', 'Y'],
    ['M', 'V', 'F', 'X', 'S', 'S', 'S', 'S', 'S', 'U', 'W', 'U', 'Z', 'P', 'W'],
    ['B', 'M', 'W', 'R', 'S', 'S', 'S', 'S', 'S', 'L', 'K', 'V', 'U', 'M', 'R'],
    ['T', 'T', 'I', 'E', 'S', 'S', 'S', 'S', 'S', 'X', 'S', 'D', 'X', 'C', 'W'],
    ['L', 'X', 'B', 'T', 'Q', 'O', 'A', 'H', 'C', 'C', 'Y', 'F', 'G', 'P', 'Z'],
    ['R', 'B', 'G', 'Q', 'G', 'V', 'U', 'W', 'O', 'E', 'P', 'I', 'O', 'Y', 'C'],
    ['Z', 'K', 'C', 'H', 'U', 'M', 'M', 'P', 'Q', 'P', 'B', 'N', 'C', 'T', 'W'],
    ['H', 'J', 'F', 'Q', 'P', 'I', 'V', 'R', 'B', 'H', 'J', 'O', 'E', 'S', 'V'],
    ['R', 'H', 'T', 'I', 'J', 'N', 'J', 'J', 'N', 'N', 'G', 'X', 'S', 'H', 'G'],
    ['L', 'C', 'I', 'G', 'B', 'L', 'Q', 'Q', 'W', 'E', 'E', 'F', 'U', 'L', 'G'],
    ['I', 'Q', 'K', 'A', 'R', 'Z', 'J', 'B', 'X', 'K', 'U', 'A', 'R', 'Z', 'V'],
    ['Q', 'W', 'J', 'Q', 'J', 'A', 'D', 'N', 'A', 'D', 'W', 'I', 'F', 'U', 'O']
]

# Find and print the bottom-right corner of the square
print(find_square(matrix))
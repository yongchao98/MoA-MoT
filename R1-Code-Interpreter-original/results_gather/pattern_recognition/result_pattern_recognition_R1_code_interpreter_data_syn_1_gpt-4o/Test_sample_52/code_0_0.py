def find_largest_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    max_size = 0
    bottom_right_corner = None

    # Iterate over each possible top-left corner of a square
    for i in range(rows):
        for j in range(cols):
            # Check possible square sizes starting from 3x3
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
                
                # If a valid square is found, update the largest square
                if is_square and size > max_size:
                    max_size = size
                    bottom_right_corner = (i + size - 1, j + size - 1)

    return bottom_right_corner

# Define the matrix
matrix = [
    ['G', 'B', 'J', 'O', 'T', 'E', 'J', 'K', 'Q', 'N', 'H', 'V'],
    ['I', 'G', 'D', 'S', 'U', 'I', 'H', 'P', 'N', 'V', 'J', 'G'],
    ['M', 'U', 'T', 'F', 'S', 'U', 'L', 'D', 'G', 'S', 'L', 'K'],
    ['S', 'B', 'F', 'N', 'R', 'L', 'F', 'Y', 'O', 'K', 'Q', 'M'],
    ['T', 'F', 'F', 'F', 'F', 'F', 'F', 'O', 'Q', 'Y', 'M', 'Z'],
    ['Z', 'F', 'F', 'F', 'F', 'F', 'F', 'O', 'W', 'U', 'E', 'C'],
    ['D', 'F', 'F', 'F', 'F', 'F', 'F', 'B', 'G', 'L', 'U', 'C'],
    ['M', 'F', 'F', 'F', 'F', 'F', 'F', 'X', 'Q', 'K', 'X', 'D'],
    ['O', 'F', 'F', 'F', 'F', 'F', 'F', 'P', 'A', 'M', 'A', 'U'],
    ['E', 'F', 'F', 'F', 'F', 'F', 'F', 'N', 'B', 'M', 'V', 'B'],
    ['N', 'N', 'U', 'X', 'U', 'G', 'U', 'S', 'S', 'N', 'Z', 'Y'],
    ['X', 'G', 'C', 'D', 'V', 'W', 'K', 'F', 'T', 'P', 'J', 'V'],
    ['M', 'F', 'K', 'Z', 'P', 'R', 'R', 'K', 'F', 'S', 'X', 'A']
]

# Find and print the bottom-right corner of the largest square
bottom_right_corner = find_largest_square(matrix)
print(bottom_right_corner)
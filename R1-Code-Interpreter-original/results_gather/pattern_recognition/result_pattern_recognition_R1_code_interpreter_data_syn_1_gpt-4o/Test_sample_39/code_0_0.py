def find_largest_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    max_size = 0
    bottom_right_corner = (-1, -1)

    # Iterate over each element to consider it as the top-left corner of a square
    for i in range(rows):
        for j in range(cols):
            # Check possible square sizes starting from 3
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
                
                # Update the largest square found
                if is_square and size > max_size:
                    max_size = size
                    bottom_right_corner = (i + size - 1, j + size - 1)

    return bottom_right_corner

# Define the matrix
matrix = [
    ['M', 'Z', 'X', 'B', 'F', 'N', 'A', 'V', 'R', 'R', 'D'],
    ['M', 'O', 'K', 'C', 'Z', 'Q', 'E', 'I', 'R', 'R', 'I'],
    ['V', 'A', 'C', 'Q', 'E', 'F', 'U', 'C', 'L', 'A', 'H'],
    ['Z', 'K', 'G', 'O', 'O', 'A', 'N', 'R', 'W', 'Y', 'Q'],
    ['K', 'Z', 'L', 'V', 'G', 'O', 'C', 'I', 'K', 'F', 'I'],
    ['G', 'E', 'F', 'E', 'C', 'O', 'Y', 'Z', 'J', 'O', 'V'],
    ['S', 'F', 'O', 'W', 'I', 'I', 'M', 'J', 'F', 'N', 'J'],
    ['G', 'V', 'T', 'M', 'M', 'M', 'M', 'M', 'B', 'I', 'X'],
    ['T', 'T', 'P', 'M', 'M', 'M', 'M', 'M', 'L', 'Q', 'Z'],
    ['E', 'K', 'Q', 'M', 'M', 'M', 'M', 'M', 'H', 'V', 'M'],
    ['D', 'I', 'L', 'M', 'M', 'M', 'M', 'M', 'B', 'T', 'G'],
    ['A', 'B', 'Q', 'M', 'M', 'M', 'M', 'M', 'N', 'W', 'F'],
    ['I', 'Q', 'R', 'K', 'G', 'T', 'Y', 'P', 'T', 'X', 'K']
]

# Find the bottom-right corner of the largest square
result = find_largest_square(matrix)
print(result)
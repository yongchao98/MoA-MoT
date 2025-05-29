def find_largest_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    max_size = 0
    bottom_right_corner = (-1, -1)

    # Start checking for squares of size 3x3 and larger
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
                
                # If a valid square is found, update the largest square info
                if is_square and size > max_size:
                    max_size = size
                    bottom_right_corner = (i + size - 1, j + size - 1)

    return bottom_right_corner

# Given matrix
matrix = [
    ['I', 'I', 'Z', 'Q', 'A', 'R', 'L', 'S', 'E', 'T', 'F', 'F', 'I', 'J', 'Z'],
    ['S', 'S', 'S', 'S', 'S', 'H', 'P', 'V', 'T', 'L', 'D', 'C', 'C', 'X', 'J'],
    ['E', 'S', 'S', 'S', 'S', 'Z', 'H', 'Q', 'G', 'H', 'Q', 'I', 'O', 'E', 'A'],
    ['R', 'S', 'S', 'S', 'S', 'E', 'Y', 'X', 'F', 'M', 'G', 'F', 'K', 'R', 'V'],
    ['X', 'S', 'S', 'S', 'S', 'A', 'K', 'O', 'M', 'Z', 'B', 'E', 'P', 'T', 'B'],
    ['A', 'Y', 'R', 'E', 'A', 'W', 'H', 'K', 'Z', 'Y', 'Q', 'E', 'Z', 'W', 'C'],
    ['N', 'F', 'B', 'E', 'C', 'K', 'H', 'M', 'W', 'A', 'L', 'R', 'X', 'P', 'K'],
    ['H', 'R', 'L', 'U', 'Q', 'B', 'G', 'I', 'U', 'E', 'Y', 'T', 'M', 'S', 'Z'],
    ['L', 'E', 'V', 'S', 'U', 'X', 'L', 'H', 'V', 'U', 'R', 'V', 'Q', 'O', 'O'],
    ['Q', 'I', 'T', 'C', 'T', 'C', 'S', 'V', 'B', 'O', 'B', 'I', 'F', 'B', 'B']
]

# Find the bottom-right corner of the largest square
result = find_largest_square(matrix)
print(result)
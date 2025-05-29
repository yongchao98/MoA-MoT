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
                
                # If a valid square is found, update the bottom-right corner
                if is_square and size > max_size:
                    max_size = size
                    bottom_right_corner = (i + size - 1, j + size - 1)

    return bottom_right_corner

# Define the matrix
matrix = [
    ['Y', 'S', 'E', 'F', 'J', 'C', 'M', 'L', 'S', 'J', 'R', 'T'],
    ['K', 'L', 'R', 'Z', 'Q', 'W', 'Q', 'K', 'H', 'F', 'X', 'N'],
    ['G', 'G', 'G', 'G', 'G', 'G', 'Q', 'P', 'S', 'O', 'F', 'X'],
    ['G', 'G', 'G', 'G', 'G', 'G', 'J', 'T', 'P', 'Z', 'Z', 'I'],
    ['G', 'G', 'G', 'G', 'G', 'G', 'L', 'M', 'Z', 'S', 'F', 'A'],
    ['G', 'G', 'G', 'G', 'G', 'G', 'Y', 'X', 'S', 'Q', 'M', 'G'],
    ['G', 'G', 'G', 'G', 'G', 'G', 'X', 'X', 'C', 'D', 'Q', 'N'],
    ['G', 'G', 'G', 'G', 'G', 'G', 'V', 'N', 'A', 'D', 'D', 'U'],
    ['J', 'L', 'Z', 'G', 'A', 'W', 'F', 'P', 'Q', 'E', 'P', 'C'],
    ['S', 'W', 'W', 'U', 'O', 'S', 'L', 'M', 'A', 'Y', 'M', 'V'],
    ['A', 'T', 'R', 'R', 'S', 'P', 'M', 'R', 'L', 'L', 'J', 'X'],
    ['N', 'M', 'K', 'T', 'W', 'H', 'Q', 'A', 'V', 'C', 'W', 'U'],
    ['I', 'C', 'V', 'P', 'R', 'L', 'Q', 'L', 'O', 'J', 'T', 'P'],
    ['M', 'Q', 'Y', 'F', 'U', 'E', 'Q', 'J', 'J', 'G', 'H', 'N'],
    ['A', 'F', 'M', 'N', 'M', 'T', 'Q', 'O', 'W', 'H', 'B', 'T']
]

# Find the bottom-right corner of the largest square
result = find_largest_square(matrix)
print(result)
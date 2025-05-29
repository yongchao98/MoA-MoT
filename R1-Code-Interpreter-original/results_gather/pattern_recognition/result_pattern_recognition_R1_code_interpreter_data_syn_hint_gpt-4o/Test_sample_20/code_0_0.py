def find_largest_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    max_side = 0
    bottom_right_corner = None

    # Iterate over each possible top-left corner of a square
    for i in range(rows):
        for j in range(cols):
            # Check possible square sizes starting from 3
            for side in range(3, min(rows - i, cols - j) + 1):
                # Check if all characters in the square are the same
                char = matrix[i][j]
                is_square = True
                for x in range(i, i + side):
                    for y in range(j, j + side):
                        if matrix[x][y] != char:
                            is_square = False
                            break
                    if not is_square:
                        break
                
                # If a valid square is found, update the largest square info
                if is_square and side > max_side:
                    max_side = side
                    bottom_right_corner = (i + side - 1, j + side - 1)

    return bottom_right_corner

# Define the matrix
matrix = [
    ['N', 'B', 'M', 'H', 'I', 'U', 'W', 'F', 'Y', 'Q', 'S', 'D', 'Z'],
    ['L', 'I', 'Z', 'T', 'C', 'E', 'R', 'W', 'H', 'A', 'Y', 'E', 'Y'],
    ['H', 'P', 'C', 'R', 'E', 'B', 'T', 'Q', 'I', 'G', 'S', 'B', 'M'],
    ['N', 'F', 'E', 'W', 'T', 'Q', 'D', 'K', 'E', 'C', 'C', 'S', 'V'],
    ['G', 'M', 'A', 'S', 'Y', 'D', 'H', 'H', 'K', 'Z', 'E', 'Q', 'V'],
    ['I', 'B', 'D', 'A', 'V', 'S', 'B', 'M', 'V', 'O', 'W', 'F', 'M'],
    ['R', 'U', 'F', 'P', 'S', 'L', 'C', 'Q', 'F', 'Y', 'P', 'V', 'T'],
    ['G', 'C', 'T', 'U', 'T', 'P', 'N', 'A', 'J', 'P', 'P', 'E', 'T'],
    ['Q', 'J', 'K', 'R', 'T', 'J', 'V', 'D', 'K', 'B', 'I', 'K', 'L'],
    ['X', 'D', 'L', 'L', 'L', 'N', 'Z', 'Q', 'J', 'U', 'M', 'K', 'S'],
    ['M', 'R', 'L', 'L', 'L', 'M', 'N', 'Y', 'S', 'X', 'F', 'U', 'Z'],
    ['F', 'H', 'L', 'L', 'L', 'K', 'R', 'N', 'Q', 'E', 'Z', 'H', 'U'],
    ['C', 'D', 'J', 'E', 'Q', 'T', 'M', 'P', 'X', 'M', 'M', 'M', 'U']
]

# Find the bottom-right corner of the largest square
result = find_largest_square(matrix)
print(result)
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
    ['T', 'E', 'H', 'Q', 'P', 'S', 'P', 'X'],
    ['Q', 'T', 'T', 'B', 'V', 'Q', 'N', 'U'],
    ['G', 'T', 'Y', 'N', 'W', 'W', 'D', 'Y'],
    ['S', 'X', 'L', 'O', 'O', 'G', 'X', 'S'],
    ['A', 'F', 'Y', 'S', 'Z', 'Z', 'K', 'A'],
    ['Z', 'N', 'T', 'A', 'D', 'K', 'J', 'M'],
    ['L', 'S', 'W', 'B', 'R', 'C', 'C', 'C'],
    ['O', 'M', 'X', 'K', 'P', 'C', 'C', 'C'],
    ['I', 'Y', 'G', 'P', 'P', 'C', 'C', 'C'],
    ['H', 'J', 'G', 'T', 'O', 'Z', 'Z', 'X'],
    ['W', 'Y', 'J', 'Z', 'G', 'H', 'N', 'D'],
    ['D', 'C', 'K', 'B', 'W', 'H', 'F', 'P'],
    ['P', 'B', 'X', 'H', 'O', 'W', 'R', 'H'],
    ['Y', 'I', 'A', 'S', 'T', 'G', 'U', 'V']
]

# Find the bottom-right corner of the largest square
result = find_largest_square(matrix)
print(result)
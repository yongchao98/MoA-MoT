def find_largest_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    max_size = 0
    bottom_right_corner = None

    # Iterate over each possible top-left corner of a square
    for i in range(rows):
        for j in range(cols):
            # Try to form a square with the top-left corner at (i, j)
            current_char = matrix[i][j]
            size = 2  # Start checking from size 3x3
            while i + size < rows and j + size < cols:
                # Check if the square of size (size+1)x(size+1) is valid
                valid_square = True
                for k in range(size + 1):
                    if matrix[i + size][j + k] != current_char or matrix[i + k][j + size] != current_char:
                        valid_square = False
                        break
                if valid_square:
                    size += 1
                else:
                    break

            # Update the largest square found
            if size > max_size:
                max_size = size
                bottom_right_corner = (i + size - 1, j + size - 1)

    return bottom_right_corner

# Define the matrix
matrix = [
    ['I', 'P', 'U', 'S', 'N', 'S', 'A'],
    ['L', 'C', 'S', 'P', 'E', 'X', 'M'],
    ['I', 'F', 'Z', 'C', 'A', 'K', 'S'],
    ['Y', 'E', 'C', 'H', 'B', 'M', 'Y'],
    ['R', 'C', 'F', 'C', 'C', 'Q', 'F'],
    ['R', 'J', 'F', 'M', 'X', 'X', 'K'],
    ['I', 'I', 'I', 'L', 'U', 'Q', 'M'],
    ['I', 'I', 'I', 'I', 'P', 'K', 'F'],
    ['I', 'I', 'I', 'M', 'U', 'P', 'D'],
    ['Q', 'L', 'F', 'A', 'R', 'W', 'K'],
    ['X', 'P', 'A', 'X', 'I', 'U', 'Z'],
    ['T', 'Z', 'U', 'X', 'P', 'B', 'F']
]

# Find the bottom-right corner of the largest square
result = find_largest_square(matrix)
print(result)
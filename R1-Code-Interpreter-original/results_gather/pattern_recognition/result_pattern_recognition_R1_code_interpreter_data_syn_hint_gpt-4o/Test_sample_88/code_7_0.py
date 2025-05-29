def find_largest_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    max_size = 0
    bottom_right_corner = (-1, -1)

    # Create a DP table to store the size of the largest square ending at each point
    dp = [[0] * cols for _ in range(rows)]

    for i in range(rows):
        for j in range(cols):
            if i == 0 or j == 0:
                # The first row and column can only have squares of size 1
                dp[i][j] = 1
            else:
                if matrix[i][j] == matrix[i-1][j] == matrix[i][j-1] == matrix[i-1][j-1]:
                    # If the current cell matches the top, left, and top-left cells
                    dp[i][j] = min(dp[i-1][j], dp[i][j-1], dp[i-1][j-1]) + 1
                else:
                    dp[i][j] = 1

            # Update the maximum size and bottom-right corner if a larger square is found
            if dp[i][j] > max_size:
                max_size = dp[i][j]
                bottom_right_corner = (i, j)

    # Adjust for the minimum size requirement of 3
    if max_size < 3:
        return None
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
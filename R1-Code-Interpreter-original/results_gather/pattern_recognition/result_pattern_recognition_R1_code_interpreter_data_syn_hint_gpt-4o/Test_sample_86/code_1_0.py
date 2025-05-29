def find_largest_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    dp = [[0] * cols for _ in range(rows)]
    max_size = 0
    bottom_right_corner = (-1, -1)

    for i in range(rows):
        for j in range(cols):
            if i == 0 or j == 0:
                dp[i][j] = 1  # The first row and column can only form squares of size 1
            else:
                if matrix[i][j] == matrix[i-1][j] == matrix[i][j-1] == matrix[i-1][j-1]:
                    dp[i][j] = min(dp[i-1][j], dp[i][j-1], dp[i-1][j-1]) + 1
                else:
                    dp[i][j] = 1

            if dp[i][j] >= 3 and dp[i][j] > max_size:
                max_size = dp[i][j]
                bottom_right_corner = (i, j)

    return bottom_right_corner

# Define the matrix
matrix = [
    ['G', 'S', 'J', 'W', 'E', 'M', 'Y', 'Y', 'C'],
    ['C', 'U', 'A', 'B', 'N', 'Z', 'Y', 'L', 'N'],
    ['E', 'L', 'F', 'F', 'F', 'F', 'Y', 'E', 'W'],
    ['Y', 'B', 'F', 'F', 'F', 'F', 'M', 'V', 'T'],
    ['H', 'P', 'F', 'F', 'F', 'F', 'R', 'N', 'H'],
    ['N', 'K', 'F', 'F', 'F', 'F', 'L', 'D', 'R'],
    ['Z', 'X', 'C', 'Y', 'I', 'W', 'D', 'V', 'W'],
    ['O', 'L', 'H', 'K', 'U', 'Y', 'K', 'U', 'P'],
    ['W', 'W', 'W', 'A', 'M', 'R', 'S', 'F', 'D'],
    ['C', 'C', 'R', 'C', 'A', 'P', 'L', 'A', 'X'],
    ['W', 'V', 'C', 'L', 'Q', 'T', 'O', 'J', 'X'],
    ['E', 'K', 'F', 'P', 'G', 'O', 'D', 'R', 'G'],
    ['G', 'A', 'I', 'B', 'I', 'A', 'M', 'A', 'C']
]

# Find and print the bottom-right corner of the largest square
print(find_largest_square(matrix))
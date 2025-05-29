def find_largest_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    dp = [[0] * cols for _ in range(rows)]
    max_size = 0
    bottom_right = (0, 0)

    for i in range(rows):
        for j in range(cols):
            if i == 0 or j == 0:
                dp[i][j] = 1
            else:
                if matrix[i][j] == matrix[i-1][j] == matrix[i][j-1] == matrix[i-1][j-1]:
                    dp[i][j] = min(dp[i-1][j], dp[i][j-1], dp[i-1][j-1]) + 1
                else:
                    dp[i][j] = 1

            if dp[i][j] > max_size:
                max_size = dp[i][j]
                bottom_right = (i, j)

    return [bottom_right[0], bottom_right[1]]

# Define the matrix
matrix = [
    ['N', 'Q', 'S', 'F', 'Z', 'W', 'T', 'N', 'M', 'S', 'H', 'G'],
    ['V', 'R', 'B', 'K', 'X', 'F', 'D', 'N', 'X', 'F', 'Q', 'S'],
    ['Z', 'B', 'N', 'U', 'L', 'I', 'V', 'C', 'V', 'G', 'W', 'B'],
    ['B', 'H', 'R', 'X', 'J', 'O', 'U', 'Q', 'E', 'O', 'T', 'G'],
    ['M', 'R', 'S', 'G', 'C', 'P', 'B', 'A', 'G', 'X', 'H', 'M'],
    ['A', 'R', 'K', 'N', 'G', 'G', 'P', 'A', 'E', 'R', 'Z', 'J'],
    ['J', 'A', 'K', 'Y', 'A', 'Z', 'I', 'V', 'E', 'R', 'B', 'X'],
    ['B', 'L', 'H', 'P', 'D', 'V', 'Y', 'X', 'N', 'B', 'G', 'K'],
    ['Q', 'U', 'B', 'G', 'S', 'Q', 'W', 'M', 'V', 'U', 'R', 'V'],
    ['B', 'C', 'L', 'L', 'H', 'X', 'K', 'B', 'U', 'H', 'T', 'N'],
    ['N', 'P', 'N', 'E', 'J', 'N', 'X', 'Z', 'J', 'S', 'L', 'K'],
    ['T', 'Q', 'J', 'D', 'T', 'W', 'N', 'K', 'Q', 'F', 'E', 'X'],
    ['K', 'X', 'I', 'C', 'O', 'L', 'S', 'M', 'M', 'M', 'B', 'Y'],
    ['X', 'J', 'W', 'Q', 'K', 'L', 'G', 'M', 'M', 'M', 'C', 'M'],
    ['V', 'D', 'T', 'N', 'X', 'L', 'I', 'M', 'M', 'M', 'T', 'A']
]

# Find and print the bottom-right corner of the largest square
print(find_largest_square(matrix))
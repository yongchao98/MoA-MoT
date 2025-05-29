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

    return bottom_right if max_size >= 3 else None

matrix = [
    ['Q', 'I', 'A', 'Y', 'W', 'B', 'D', 'K', 'K', 'E', 'J', 'F'],
    ['Q', 'R', 'M', 'U', 'E', 'Q', 'I', 'I', 'K', 'S', 'X', 'K'],
    ['Y', 'D', 'R', 'E', 'F', 'U', 'I', 'W', 'H', 'V', 'Y', 'I'],
    ['L', 'Z', 'S', 'E', 'F', 'P', 'A', 'M', 'D', 'M', 'Y', 'U'],
    ['X', 'P', 'C', 'C', 'C', 'C', 'H', 'Z', 'V', 'O', 'V', 'V'],
    ['L', 'T', 'C', 'C', 'C', 'C', 'P', 'Y', 'P', 'A', 'Q', 'V'],
    ['X', 'O', 'C', 'C', 'C', 'C', 'P', 'C', 'S', 'X', 'O', 'B'],
    ['P', 'N', 'C', 'C', 'C', 'C', 'Z', 'E', 'A', 'X', 'B', 'L'],
    ['K', 'R', 'U', 'A', 'O', 'Z', 'C', 'I', 'G', 'S', 'L', 'E'],
    ['D', 'G', 'O', 'Q', 'O', 'S', 'I', 'U', 'P', 'Z', 'E', 'N']
]

bottom_right_corner = find_largest_square(matrix)
print(bottom_right_corner)
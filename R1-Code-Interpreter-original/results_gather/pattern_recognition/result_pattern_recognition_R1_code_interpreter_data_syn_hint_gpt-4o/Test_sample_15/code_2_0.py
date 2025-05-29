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
                    # If the current cell matches the top, left, and top-left cells, extend the square
                    dp[i][j] = min(dp[i-1][j], dp[i][j-1], dp[i-1][j-1]) + 1
                else:
                    dp[i][j] = 1

            # Update the maximum size and bottom-right corner
            if dp[i][j] > max_size and dp[i][j] >= 3:
                max_size = dp[i][j]
                bottom_right_corner = (i, j)

    return bottom_right_corner

matrix = [
    ['S', 'P', 'P', 'P', 'P', 'W', 'V', 'L', 'Z', 'N', 'W', 'Y'],
    ['P', 'P', 'P', 'P', 'P', 'Q', 'G', 'C', 'F', 'C', 'S', 'U'],
    ['G', 'P', 'P', 'P', 'P', 'F', 'D', 'K', 'V', 'W', 'W', 'Y'],
    ['A', 'P', 'P', 'P', 'P', 'U', 'F', 'C', 'I', 'D', 'W', 'Q'],
    ['A', 'Z', 'E', 'B', 'Z', 'A', 'W', 'A', 'E', 'W', 'C', 'E'],
    ['C', 'V', 'L', 'X', 'R', 'G', 'W', 'H', 'J', 'N', 'T', 'J'],
    ['F', 'F', 'H', 'R', 'P', 'A', 'J', 'F', 'L', 'S', 'U', 'W'],
    ['Z', 'T', 'C', 'V', 'D', 'U', 'J', 'E', 'A', 'V', 'B', 'B'],
    ['A', 'L', 'T', 'Z', 'Q', 'F', 'S', 'J', 'A', 'J', 'G', 'J']
]

result = find_largest_square(matrix)
print(result)
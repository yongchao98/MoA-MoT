def find_largest_square(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    max_size = 0
    bottom_right_corner = None

    # Create a DP table to store the size of the largest square ending at each position
    dp = [[0] * cols for _ in range(rows)]

    for i in range(rows):
        for j in range(cols):
            if i == 0 or j == 0:
                # The first row and column can only have squares of size 1
                dp[i][j] = 1
            else:
                if matrix[i][j] == matrix[i-1][j] == matrix[i][j-1] == matrix[i-1][j-1]:
                    dp[i][j] = min(dp[i-1][j], dp[i][j-1], dp[i-1][j-1]) + 1
                else:
                    dp[i][j] = 1

            # Update the maximum size and bottom-right corner
            if dp[i][j] > max_size:
                max_size = dp[i][j]
                bottom_right_corner = [i, j]

    # Adjust for the minimum size requirement of 3
    if max_size < 3:
        return None

    return bottom_right_corner

matrix = [
    ['L', 'H', 'P', 'T', 'T', 'K', 'A', 'S'],
    ['B', 'W', 'W', 'W', 'D', 'I', 'J', 'W'],
    ['U', 'W', 'W', 'W', 'A', 'A', 'E', 'A'],
    ['O', 'W', 'W', 'W', 'T', 'T', 'Y', 'S'],
    ['D', 'M', 'Z', 'V', 'D', 'Z', 'Y', 'V'],
    ['S', 'C', 'T', 'W', 'H', 'C', 'Z', 'P'],
    ['K', 'M', 'N', 'Q', 'E', 'H', 'D', 'N'],
    ['E', 'K', 'S', 'X', 'A', 'U', 'I', 'D'],
    ['F', 'F', 'D', 'S', 'V', 'J', 'P', 'J'],
    ['I', 'I', 'S', 'I', 'T', 'I', 'V', 'Z'],
    ['R', 'Z', 'L', 'W', 'Y', 'J', 'D', 'E']
]

result = find_largest_square(matrix)
print(result)
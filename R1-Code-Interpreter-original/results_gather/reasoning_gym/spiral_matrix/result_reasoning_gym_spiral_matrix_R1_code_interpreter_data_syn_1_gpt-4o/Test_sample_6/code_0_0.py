def spiral_order(matrix):
    if not matrix:
        return []

    top, bottom = 0, len(matrix) - 1
    left, right = 0, len(matrix[0]) - 1
    result = []

    while top <= bottom and left <= right:
        # Traverse from left to right across the top row
        for i in range(left, right + 1):
            result.append(matrix[top][i])
        top += 1

        # Traverse from top to bottom down the right column
        for i in range(top, bottom + 1):
            result.append(matrix[i][right])
        right -= 1

        if top <= bottom:
            # Traverse from right to left across the bottom row
            for i in range(right, left - 1, -1):
                result.append(matrix[bottom][i])
            bottom -= 1

        if left <= right:
            # Traverse from bottom to top up the left column
            for i in range(bottom, top - 1, -1):
                result.append(matrix[i][left])
            left += 1

    return result

matrix = [
    [9, 5, 1, 5, 9, 9, 9, 8],
    [6, 4, 5, 4, 3, 6, 5, 6],
    [7, 2, 6, 1, 4, 0, 2, 4],
    [1, 6, 8, 3, 7, 2, 2, 1],
    [1, 6, 7, 0, 4, 2, 4, 5],
    [8, 0, 3, 2, 9, 1, 1, 8],
    [7, 4, 8, 5, 1, 7, 2, 7],
    [7, 4, 5, 3, 9, 6, 9, 5]
]

result = spiral_order(matrix)
print("<<<" + " ".join(map(str, result)) + ">>>")
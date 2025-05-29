def spiral_order(matrix):
    if not matrix:
        return []

    result = []
    top, bottom = 0, len(matrix) - 1
    left, right = 0, len(matrix[0]) - 1

    while top <= bottom and left <= right:
        # Traverse from left to right
        for i in range(left, right + 1):
            result.append(matrix[top][i])
        top += 1

        # Traverse from top to bottom
        for i in range(top, bottom + 1):
            result.append(matrix[i][right])
        right -= 1

        if top <= bottom:
            # Traverse from right to left
            for i in range(right, left - 1, -1):
                result.append(matrix[bottom][i])
            bottom -= 1

        if left <= right:
            # Traverse from bottom to top
            for i in range(bottom, top - 1, -1):
                result.append(matrix[i][left])
            left += 1

    return result

matrix = [
    [4, 5, 9, 2, 2, 1, 0, 8, 8, 1],
    [5, 4, 9, 6, 6, 2, 0, 1, 7, 5],
    [6, 9, 8, 9, 5, 9, 5, 3, 1, 5],
    [6, 3, 4, 4, 6, 5, 8, 6, 5, 6],
    [3, 7, 2, 8, 6, 8, 9, 1, 1, 4],
    [7, 9, 0, 4, 0, 4, 9, 4, 2, 2],
    [4, 8, 9, 7, 3, 6, 4, 1, 2, 4],
    [2, 1, 1, 3, 0, 6, 0, 4, 9, 2],
    [4, 5, 7, 7, 1, 5, 2, 1, 2, 6],
    [6, 8, 8, 1, 5, 3, 3, 3, 7, 6]
]

result = spiral_order(matrix)
print("<<<" + " ".join(map(str, result)) + ">>>")
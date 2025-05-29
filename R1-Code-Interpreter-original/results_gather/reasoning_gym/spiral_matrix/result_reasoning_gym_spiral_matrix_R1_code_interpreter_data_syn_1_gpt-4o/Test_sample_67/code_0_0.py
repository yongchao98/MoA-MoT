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
    [8, 5, 3, 4, 7, 2, 5, 9],
    [2, 1, 2, 9, 7, 7, 2, 8],
    [9, 3, 8, 6, 0, 5, 3, 6],
    [9, 3, 3, 9, 8, 7, 7, 8],
    [7, 2, 1, 9, 0, 8, 0, 1],
    [3, 2, 2, 0, 6, 1, 9, 0],
    [2, 1, 2, 6, 4, 6, 7, 6],
    [3, 3, 6, 3, 5, 5, 1, 8]
]

spiral_result = spiral_order(matrix)
print("<<<" + " ".join(map(str, spiral_result)) + ">>>")
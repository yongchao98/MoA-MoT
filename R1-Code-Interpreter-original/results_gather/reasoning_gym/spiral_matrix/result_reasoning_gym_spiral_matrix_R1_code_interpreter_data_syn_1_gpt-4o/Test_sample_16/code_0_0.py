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

        # Traverse downwards
        for i in range(top, bottom + 1):
            result.append(matrix[i][right])
        right -= 1

        if top <= bottom:
            # Traverse from right to left
            for i in range(right, left - 1, -1):
                result.append(matrix[bottom][i])
            bottom -= 1

        if left <= right:
            # Traverse upwards
            for i in range(bottom, top - 1, -1):
                result.append(matrix[i][left])
            left += 1

    return result

matrix = [
    [0, 2, 7, 2, 8, 5, 1, 9, 4, 6],
    [9, 9, 5, 5, 0, 1, 7, 1, 3, 4],
    [2, 4, 4, 4, 2, 6, 3, 3, 5, 9],
    [1, 9, 2, 3, 4, 6, 2, 4, 2, 0],
    [7, 8, 8, 2, 3, 0, 3, 2, 9, 0],
    [8, 5, 6, 5, 2, 8, 8, 3, 6, 6],
    [4, 6, 3, 2, 8, 3, 5, 9, 6, 1],
    [4, 6, 5, 7, 2, 2, 0, 7, 6, 7],
    [9, 6, 1, 1, 1, 9, 9, 7, 4, 6],
    [9, 2, 5, 8, 8, 3, 7, 8, 0, 5]
]

spiral_result = spiral_order(matrix)
print("<<<" + " ".join(map(str, spiral_result)) + ">>>")
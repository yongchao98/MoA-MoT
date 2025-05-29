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
    [6, 7, 3, 9, 3, 7, 3, 4, 6],
    [7, 1, 5, 6, 6, 9, 1, 6, 7],
    [3, 2, 1, 0, 5, 5, 6, 5, 6],
    [3, 2, 2, 8, 8, 6, 6, 4, 0],
    [5, 9, 6, 0, 9, 9, 1, 4, 8],
    [8, 3, 0, 8, 2, 8, 9, 7, 9],
    [1, 1, 7, 8, 8, 7, 2, 4, 1],
    [0, 6, 2, 0, 5, 0, 1, 9, 6],
    [7, 0, 8, 7, 8, 6, 2, 4, 9]
]

spiral = spiral_order(matrix)
print("<<<" + " ".join(map(str, spiral)) + ">>>")
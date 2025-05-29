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
    [2, 8, 1, 6, 1, 2, 9, 7, 2, 4],
    [8, 3, 0, 1, 7, 2, 5, 0, 3, 9],
    [5, 0, 5, 8, 4, 8, 6, 5, 6, 0],
    [0, 8, 2, 0, 8, 7, 7, 3, 2, 6],
    [4, 0, 8, 5, 1, 1, 9, 6, 3, 8],
    [6, 0, 8, 2, 1, 8, 5, 9, 5, 9],
    [7, 7, 5, 6, 3, 7, 6, 7, 5, 4],
    [8, 2, 8, 3, 1, 2, 2, 6, 4, 8],
    [2, 6, 4, 8, 9, 5, 0, 5, 4, 5],
    [4, 2, 1, 6, 6, 5, 1, 9, 7, 6]
]

spiral = spiral_order(matrix)
print("<<<" + " ".join(map(str, spiral)) + ">>>")
def spiral_order(matrix):
    if not matrix:
        return []

    result = []
    top, bottom = 0, len(matrix) - 1
    left, right = 0, len(matrix[0]) - 1

    while top <= bottom and left <= right:
        # Traverse from left to right along the top boundary
        for i in range(left, right + 1):
            result.append(matrix[top][i])
        top += 1

        # Traverse from top to bottom along the right boundary
        for i in range(top, bottom + 1):
            result.append(matrix[i][right])
        right -= 1

        if top <= bottom:
            # Traverse from right to left along the bottom boundary
            for i in range(right, left - 1, -1):
                result.append(matrix[bottom][i])
            bottom -= 1

        if left <= right:
            # Traverse from bottom to top along the left boundary
            for i in range(bottom, top - 1, -1):
                result.append(matrix[i][left])
            left += 1

    return result

matrix = [
    [1, 9, 9, 5, 2, 9, 7, 3],
    [1, 1, 5, 0, 7, 0, 4, 9],
    [3, 5, 4, 7, 8, 4, 3, 4],
    [6, 5, 3, 3, 2, 7, 1, 9],
    [6, 7, 7, 0, 1, 4, 1, 8],
    [8, 2, 5, 9, 0, 1, 4, 0],
    [2, 1, 5, 5, 6, 4, 0, 3],
    [1, 6, 6, 0, 2, 8, 8, 5]
]

spiral_result = spiral_order(matrix)
print("<<<" + " ".join(map(str, spiral_result)) + ">>>")
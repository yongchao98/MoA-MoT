def spiral_order(matrix):
    if not matrix:
        return []

    result = []
    top, bottom = 0, len(matrix) - 1
    left, right = 0, len(matrix[0]) - 1

    while top <= bottom and left <= right:
        # Traverse from left to right along the top row
        for i in range(left, right + 1):
            result.append(matrix[top][i])
        top += 1

        # Traverse downwards along the right column
        for i in range(top, bottom + 1):
            result.append(matrix[i][right])
        right -= 1

        if top <= bottom:
            # Traverse from right to left along the bottom row
            for i in range(right, left - 1, -1):
                result.append(matrix[bottom][i])
            bottom -= 1

        if left <= right:
            # Traverse upwards along the left column
            for i in range(bottom, top - 1, -1):
                result.append(matrix[i][left])
            left += 1

    return result

matrix = [
    [7, 0, 8, 5, 2, 7, 0, 9],
    [2, 3, 6, 0, 0, 5, 2, 3],
    [7, 5, 3, 2, 5, 8, 7, 5],
    [6, 8, 3, 2, 2, 2, 6, 3],
    [6, 1, 5, 8, 4, 7, 7, 2],
    [6, 6, 2, 2, 3, 1, 6, 8],
    [7, 1, 4, 4, 1, 0, 8, 8],
    [6, 2, 0, 7, 2, 3, 1, 0]
]

spiral = spiral_order(matrix)
print("<<<" + " ".join(map(str, spiral)) + ">>>")
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
    [1, 8, 8, 1, 2, 3, 8, 6, 0, 8],
    [4, 7, 6, 1, 7, 5, 7, 0, 6, 9],
    [0, 6, 9, 5, 6, 6, 4, 2, 8, 7],
    [4, 5, 6, 1, 8, 6, 7, 7, 2, 0],
    [7, 0, 0, 8, 4, 8, 6, 7, 5, 0],
    [4, 4, 3, 3, 0, 5, 7, 2, 2, 2],
    [0, 6, 1, 9, 9, 9, 4, 4, 5, 3],
    [3, 4, 6, 8, 7, 0, 7, 5, 0, 7],
    [8, 5, 8, 4, 2, 6, 9, 9, 5, 3],
    [1, 5, 1, 8, 6, 9, 7, 0, 1, 6]
]

spiral_result = spiral_order(matrix)
print("<<<" + " ".join(map(str, spiral_result)) + ">>>")
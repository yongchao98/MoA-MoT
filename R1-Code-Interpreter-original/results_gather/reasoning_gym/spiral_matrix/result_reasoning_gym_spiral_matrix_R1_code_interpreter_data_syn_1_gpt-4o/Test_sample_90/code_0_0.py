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
    [6, 7, 1, 9, 6, 1, 7, 6, 8, 6],
    [1, 0, 5, 8, 1, 2, 4, 3, 5, 7],
    [9, 9, 6, 5, 7, 5, 3, 3, 2, 1],
    [2, 9, 1, 0, 9, 0, 7, 6, 7, 7],
    [6, 3, 8, 4, 4, 6, 5, 5, 1, 2],
    [1, 0, 5, 5, 9, 8, 4, 0, 2, 1],
    [5, 8, 2, 3, 3, 7, 9, 7, 6, 3],
    [9, 4, 4, 9, 5, 6, 2, 3, 5, 1],
    [7, 2, 1, 9, 5, 9, 9, 7, 3, 6],
    [3, 0, 7, 1, 9, 6, 4, 0, 8, 7]
]

spiral = spiral_order(matrix)
print("<<<" + " ".join(map(str, spiral)) + ">>>")
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
    [3, 9, 2, 0, 5, 9, 3, 0, 4, 2],
    [8, 7, 9, 9, 2, 8, 7, 2, 5, 9],
    [3, 5, 8, 9, 7, 3, 7, 0, 6, 2],
    [3, 1, 1, 8, 7, 0, 5, 1, 9, 3],
    [9, 6, 1, 6, 4, 2, 0, 9, 7, 6],
    [3, 3, 4, 2, 3, 1, 4, 3, 9, 4],
    [6, 5, 1, 8, 0, 3, 2, 5, 5, 1],
    [5, 4, 0, 9, 8, 8, 8, 6, 9, 2],
    [9, 3, 6, 0, 6, 0, 8, 2, 5, 1],
    [7, 1, 2, 8, 6, 9, 2, 4, 4, 2]
]

spiral = spiral_order(matrix)
print("<<<" + " ".join(map(str, spiral)) + ">>>")
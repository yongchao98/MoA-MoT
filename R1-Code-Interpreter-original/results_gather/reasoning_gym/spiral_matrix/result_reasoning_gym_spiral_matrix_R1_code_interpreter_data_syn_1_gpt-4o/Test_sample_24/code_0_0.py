def spiral_order(matrix):
    if not matrix or not matrix[0]:
        return []

    top, bottom = 0, len(matrix) - 1
    left, right = 0, len(matrix[0]) - 1
    result = []

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
    [0, 8, 5, 6, 7, 0, 9, 4, 3, 7],
    [2, 0, 1, 5, 0, 6, 8, 5, 2, 7],
    [7, 2, 8, 8, 8, 2, 5, 7, 6, 9],
    [1, 1, 8, 4, 0, 6, 2, 2, 2, 9],
    [5, 4, 2, 8, 8, 2, 0, 8, 0, 3],
    [8, 1, 1, 0, 9, 4, 3, 7, 6, 0],
    [0, 3, 7, 7, 9, 3, 2, 5, 3, 9],
    [0, 3, 9, 3, 3, 9, 7, 9, 7, 5],
    [8, 6, 8, 0, 8, 0, 2, 3, 4, 9],
    [2, 3, 9, 4, 8, 0, 1, 2, 9, 3]
]

spiral = spiral_order(matrix)
print("<<<" + " ".join(map(str, spiral)) + ">>>")
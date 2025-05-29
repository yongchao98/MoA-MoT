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
    [3, 6, 3, 2, 3, 0, 4, 0, 7],
    [5, 4, 0, 8, 1, 4, 4, 9, 0],
    [1, 0, 8, 6, 6, 2, 5, 3, 3],
    [2, 3, 4, 5, 9, 8, 9, 8, 6],
    [1, 5, 5, 1, 4, 5, 3, 9, 3],
    [1, 6, 2, 4, 7, 0, 2, 8, 6],
    [2, 9, 9, 0, 0, 7, 7, 4, 5],
    [6, 6, 5, 9, 2, 1, 5, 0, 2],
    [5, 0, 9, 1, 0, 1, 7, 7, 1]
]

spiral = spiral_order(matrix)
print("<<<" + " ".join(map(str, spiral)) + ">>>")
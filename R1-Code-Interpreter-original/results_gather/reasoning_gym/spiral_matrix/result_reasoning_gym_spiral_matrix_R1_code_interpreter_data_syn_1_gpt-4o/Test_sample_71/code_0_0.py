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
    [5, 0, 6, 3, 8, 4, 3, 6, 1],
    [0, 0, 1, 7, 1, 7, 9, 9, 7],
    [1, 4, 4, 0, 1, 6, 7, 1, 7],
    [1, 1, 9, 8, 0, 9, 9, 5, 6],
    [1, 4, 2, 7, 6, 3, 7, 8, 3],
    [8, 9, 4, 3, 9, 1, 4, 8, 3],
    [4, 2, 5, 9, 0, 7, 6, 3, 2],
    [7, 6, 9, 3, 8, 4, 9, 1, 2],
    [2, 8, 2, 1, 3, 4, 7, 5, 9]
]

spiral = spiral_order(matrix)
print("<<<" + " ".join(map(str, spiral)) + ">>>")
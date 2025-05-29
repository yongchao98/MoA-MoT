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
    [2, 4, 4, 3, 7, 4, 6, 2, 9, 6],
    [2, 0, 3, 0, 6, 3, 4, 7, 6, 8],
    [0, 9, 8, 7, 5, 5, 1, 5, 8, 5],
    [4, 6, 0, 0, 0, 8, 6, 3, 5, 2],
    [8, 9, 5, 3, 2, 3, 6, 5, 0, 4],
    [2, 9, 6, 5, 6, 1, 9, 5, 8, 1],
    [1, 2, 8, 0, 4, 6, 0, 1, 6, 3],
    [0, 4, 4, 4, 3, 7, 7, 8, 4, 3],
    [2, 1, 5, 2, 3, 3, 7, 8, 9, 6],
    [7, 9, 5, 9, 7, 0, 3, 2, 7, 1]
]

spiral_result = spiral_order(matrix)
print("<<<" + " ".join(map(str, spiral_result)) + ">>>")
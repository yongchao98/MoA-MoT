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
    [7, 7, 3, 6, 4, 3, 9],
    [9, 0, 3, 3, 4, 8, 5],
    [2, 1, 0, 8, 5, 6, 1],
    [9, 5, 3, 0, 4, 7, 6],
    [0, 4, 8, 4, 6, 2, 7],
    [9, 6, 1, 1, 1, 3, 6],
    [0, 9, 7, 3, 7, 1, 0]
]

result = spiral_order(matrix)
print("<<<" + " ".join(map(str, result)) + ">>>")
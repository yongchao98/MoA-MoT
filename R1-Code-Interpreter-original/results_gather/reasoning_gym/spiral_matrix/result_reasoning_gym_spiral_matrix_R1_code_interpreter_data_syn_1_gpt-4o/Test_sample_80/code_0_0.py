def spiral_order(matrix):
    if not matrix:
        return []

    result = []
    top, bottom = 0, len(matrix) - 1
    left, right = 0, len(matrix[0]) - 1

    while top <= bottom and left <= right:
        # Traverse from left to right on the top row
        for i in range(left, right + 1):
            result.append(matrix[top][i])
        top += 1

        # Traverse from top to bottom on the right column
        for i in range(top, bottom + 1):
            result.append(matrix[i][right])
        right -= 1

        if top <= bottom:
            # Traverse from right to left on the bottom row
            for i in range(right, left - 1, -1):
                result.append(matrix[bottom][i])
            bottom -= 1

        if left <= right:
            # Traverse from bottom to top on the left column
            for i in range(bottom, top - 1, -1):
                result.append(matrix[i][left])
            left += 1

    return result

matrix = [
    [9, 1, 6, 8, 9, 4, 5, 7, 8],
    [3, 8, 3, 0, 3, 5, 9, 3, 2],
    [9, 1, 7, 0, 4, 2, 4, 1, 0],
    [4, 9, 5, 2, 9, 1, 7, 9, 3],
    [6, 8, 4, 4, 7, 6, 2, 7, 7],
    [8, 9, 1, 9, 3, 6, 1, 1, 6],
    [8, 7, 6, 4, 6, 1, 8, 3, 8],
    [5, 6, 0, 3, 9, 3, 5, 9, 6],
    [1, 2, 7, 1, 2, 9, 8, 4, 9]
]

spiral = spiral_order(matrix)
print("<<<" + " ".join(map(str, spiral)) + ">>>")
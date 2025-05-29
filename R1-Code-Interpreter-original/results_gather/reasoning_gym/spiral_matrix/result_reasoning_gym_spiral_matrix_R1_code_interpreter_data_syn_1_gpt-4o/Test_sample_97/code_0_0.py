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

        # Traverse from top to bottom along the right column
        for i in range(top, bottom + 1):
            result.append(matrix[i][right])
        right -= 1

        if top <= bottom:
            # Traverse from right to left along the bottom row
            for i in range(right, left - 1, -1):
                result.append(matrix[bottom][i])
            bottom -= 1

        if left <= right:
            # Traverse from bottom to top along the left column
            for i in range(bottom, top - 1, -1):
                result.append(matrix[i][left])
            left += 1

    return result

matrix = [
    [4, 6, 3, 3, 8, 6, 5, 1, 3, 9],
    [9, 4, 5, 8, 2, 3, 3, 8, 9, 2],
    [5, 1, 1, 7, 7, 0, 4, 4, 7, 3],
    [0, 2, 3, 9, 4, 6, 3, 3, 9, 8],
    [9, 2, 1, 0, 2, 6, 3, 2, 1, 4],
    [2, 7, 1, 2, 0, 2, 5, 3, 9, 2],
    [6, 6, 9, 2, 5, 4, 3, 2, 8, 8],
    [5, 8, 0, 6, 2, 1, 4, 4, 9, 5],
    [6, 8, 0, 8, 3, 2, 9, 0, 7, 7],
    [9, 2, 5, 0, 4, 4, 6, 2, 9, 7]
]

spiral = spiral_order(matrix)
print("<<<" + " ".join(map(str, spiral)) + ">>>")
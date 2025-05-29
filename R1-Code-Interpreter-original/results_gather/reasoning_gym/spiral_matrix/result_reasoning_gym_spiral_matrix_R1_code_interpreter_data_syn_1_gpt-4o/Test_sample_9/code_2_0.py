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
    [4, 0, 1, 1, 5, 2, 0, 1],
    [4, 9, 6, 6, 7, 1, 5, 9],
    [4, 0, 0, 2, 4, 2, 0, 0],
    [7, 1, 3, 6, 7, 8, 0, 1],
    [9, 1, 2, 2, 9, 6, 0, 0],
    [4, 4, 6, 9, 8, 0, 2, 2],
    [9, 9, 6, 3, 9, 5, 3, 5],
    [7, 3, 2, 2, 6, 1, 3, 2]
]

spiral_result = spiral_order(matrix)
print("<<<" + " ".join(map(str, spiral_result)) + ">>>")
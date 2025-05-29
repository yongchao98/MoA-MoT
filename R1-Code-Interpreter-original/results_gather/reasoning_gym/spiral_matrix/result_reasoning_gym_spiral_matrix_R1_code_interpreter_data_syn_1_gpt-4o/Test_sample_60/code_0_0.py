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
    [7, 4, 8, 2, 1, 4, 0, 2],
    [9, 4, 6, 0, 7, 9, 7, 7],
    [1, 2, 7, 2, 6, 1, 3, 9],
    [1, 0, 3, 2, 2, 9, 8, 0],
    [9, 2, 4, 9, 6, 9, 4, 1],
    [5, 2, 9, 6, 7, 5, 9, 8],
    [8, 1, 1, 5, 0, 5, 2, 2],
    [0, 6, 4, 8, 9, 9, 9, 2]
]

spiral = spiral_order(matrix)
print("<<<" + " ".join(map(str, spiral)) + ">>>")
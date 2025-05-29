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
    [4, 3, 2, 9, 3, 0, 4, 7, 3, 8],
    [7, 5, 5, 1, 5, 1, 4, 6, 9, 9],
    [8, 7, 6, 7, 2, 4, 2, 2, 1, 9],
    [4, 2, 1, 0, 3, 1, 2, 6, 2, 4],
    [3, 1, 0, 6, 0, 7, 1, 7, 8, 2],
    [0, 2, 5, 7, 2, 5, 6, 9, 5, 3],
    [9, 7, 9, 4, 2, 9, 9, 0, 0, 0],
    [3, 7, 6, 8, 0, 8, 4, 3, 3, 0],
    [1, 4, 8, 2, 2, 7, 5, 7, 2, 3],
    [2, 8, 6, 8, 7, 5, 4, 9, 4, 8]
]

spiral = spiral_order(matrix)
print("<<<" + " ".join(map(str, spiral)) + ">>>")
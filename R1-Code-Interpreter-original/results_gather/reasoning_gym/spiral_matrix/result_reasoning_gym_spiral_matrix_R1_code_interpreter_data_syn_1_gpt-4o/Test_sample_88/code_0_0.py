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
    [0, 2, 4, 9, 7, 9, 0, 8, 9, 9],
    [6, 3, 1, 9, 6, 5, 2, 1, 3, 8],
    [8, 3, 0, 1, 4, 3, 5, 7, 1, 3],
    [5, 4, 6, 6, 1, 5, 2, 4, 5, 0],
    [5, 6, 5, 8, 3, 4, 8, 6, 5, 6],
    [3, 4, 6, 4, 8, 2, 0, 7, 2, 7],
    [5, 5, 3, 4, 6, 9, 8, 1, 6, 5],
    [0, 8, 3, 7, 4, 0, 7, 2, 4, 3],
    [1, 6, 6, 3, 8, 0, 5, 9, 1, 5],
    [1, 9, 4, 5, 3, 2, 8, 3, 3, 3]
]

spiral = spiral_order(matrix)
print("<<<" + " ".join(map(str, spiral)) + ">>>")
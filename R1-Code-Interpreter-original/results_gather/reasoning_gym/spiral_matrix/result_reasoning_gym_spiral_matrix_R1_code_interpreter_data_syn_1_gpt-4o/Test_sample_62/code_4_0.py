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
    [3, 5, 0, 6, 3, 9, 0, 8, 8, 6],
    [6, 2, 8, 2, 0, 2, 5, 4, 5, 2],
    [0, 6, 4, 6, 2, 7, 3, 7, 5, 0],
    [0, 0, 5, 2, 9, 6, 0, 5, 2, 0],
    [4, 8, 2, 0, 8, 9, 4, 9, 4, 8],
    [4, 5, 8, 4, 7, 7, 4, 3, 0, 9],
    [5, 9, 7, 3, 9, 2, 5, 8, 0, 7],
    [9, 7, 1, 4, 2, 5, 2, 8, 7, 3],
    [1, 4, 1, 9, 2, 3, 6, 7, 3, 6],
    [5, 4, 5, 2, 3, 1, 6, 6, 2, 3]
]

spiral = spiral_order(matrix)
print("<<<" + " ".join(map(str, spiral)) + ">>>")
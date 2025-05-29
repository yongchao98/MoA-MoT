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
    [4, 2, 4, 8, 0, 1, 4, 7, 9, 5],
    [5, 1, 8, 5, 6, 0, 0, 0, 3, 6],
    [3, 9, 5, 1, 8, 0, 7, 7, 6, 8],
    [5, 8, 6, 3, 1, 6, 4, 5, 8, 2],
    [0, 2, 8, 0, 1, 7, 8, 6, 3, 1],
    [5, 9, 9, 0, 1, 1, 8, 4, 7, 0],
    [9, 2, 5, 2, 7, 8, 7, 9, 8, 0],
    [0, 9, 2, 0, 2, 5, 4, 9, 2, 1],
    [7, 2, 5, 5, 6, 8, 7, 9, 1, 4],
    [9, 7, 9, 2, 6, 4, 1, 1, 7, 0]
]

spiral = spiral_order(matrix)
print("<<<" + " ".join(map(str, spiral)) + ">>>")
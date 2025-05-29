def spiral_order(matrix):
    if not matrix:
        return []

    result = []
    top, bottom = 0, len(matrix) - 1
    left, right = 0, len(matrix[0]) - 1

    while top <= bottom and left <= right:
        # Traverse from left to right across the top row
        for i in range(left, right + 1):
            result.append(matrix[top][i])
        top += 1

        # Traverse down the rightmost column
        for i in range(top, bottom + 1):
            result.append(matrix[i][right])
        right -= 1

        if top <= bottom:
            # Traverse from right to left across the bottom row
            for i in range(right, left - 1, -1):
                result.append(matrix[bottom][i])
            bottom -= 1

        if left <= right:
            # Traverse up the leftmost column
            for i in range(bottom, top - 1, -1):
                result.append(matrix[i][left])
            left += 1

    return result

matrix = [
    [5, 9, 9, 9, 2, 3, 5, 3, 0],
    [8, 5, 9, 7, 1, 1, 2, 9, 8],
    [1, 9, 6, 5, 3, 1, 5, 2, 2],
    [1, 7, 6, 3, 8, 7, 7, 9, 5],
    [8, 7, 0, 9, 4, 1, 8, 5, 5],
    [6, 0, 8, 1, 4, 3, 8, 9, 3],
    [6, 2, 9, 5, 4, 4, 4, 3, 0],
    [3, 3, 8, 6, 5, 6, 2, 5, 7],
    [4, 1, 1, 8, 6, 1, 5, 1, 8]
]

spiral_result = spiral_order(matrix)
print("<<<" + " ".join(map(str, spiral_result)) + ">>>")
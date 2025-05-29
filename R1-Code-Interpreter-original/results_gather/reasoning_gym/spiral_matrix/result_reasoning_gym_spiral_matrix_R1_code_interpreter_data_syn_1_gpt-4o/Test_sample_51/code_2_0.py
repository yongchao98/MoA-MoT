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
    [6, 7, 4, 8, 7, 0, 7, 2, 6],
    [8, 2, 5, 2, 2, 6, 9, 1, 4],
    [7, 1, 2, 0, 2, 3, 7, 5, 5],
    [1, 7, 0, 9, 9, 9, 5, 8, 7],
    [0, 7, 9, 5, 0, 6, 1, 4, 1],
    [3, 7, 9, 0, 9, 2, 3, 4, 6],
    [6, 6, 9, 7, 7, 0, 0, 8, 3],
    [1, 3, 7, 4, 5, 3, 3, 2, 9],
    [2, 8, 2, 2, 2, 5, 1, 5, 4]
]

spiral = spiral_order(matrix)
print("<<<" + " ".join(map(str, spiral)) + ">>>")
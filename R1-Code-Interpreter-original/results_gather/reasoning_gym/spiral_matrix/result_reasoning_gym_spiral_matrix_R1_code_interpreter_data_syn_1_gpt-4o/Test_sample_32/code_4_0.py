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
    [5, 8, 4, 3, 8, 0, 6, 4, 8],
    [3, 3, 7, 1, 1, 4, 6, 3, 9],
    [1, 9, 6, 4, 6, 2, 7, 8, 3],
    [0, 8, 2, 3, 0, 2, 1, 2, 0],
    [5, 1, 6, 1, 2, 7, 4, 6, 3],
    [7, 2, 3, 0, 0, 1, 3, 9, 3],
    [4, 3, 3, 8, 4, 2, 2, 5, 1],
    [1, 9, 0, 2, 0, 4, 0, 4, 7],
    [2, 3, 9, 0, 4, 4, 2, 0, 3]
]

spiral_result = spiral_order(matrix)
print("<<<" + " ".join(map(str, spiral_result)) + ">>>")
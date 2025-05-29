def spiral_order(matrix):
    if not matrix:
        return []

    result = []
    while matrix:
        # Take the first row
        result += matrix.pop(0)
        # Take the last element of each remaining row
        if matrix and matrix[0]:
            for row in matrix:
                result.append(row.pop())
        # Take the last row in reverse order
        if matrix:
            result += matrix.pop()[::-1]
        # Take the first element of each remaining row (in reverse order)
        if matrix and matrix[0]:
            for row in matrix[::-1]:
                result.append(row.pop(0))
    return result

matrix = [
    [9, 2, 1, 7, 0, 5, 3, 2, 4],
    [9, 1, 8, 3, 6, 1, 0, 5, 4],
    [2, 0, 0, 4, 8, 6, 0, 9, 9],
    [4, 9, 1, 8, 8, 6, 0, 1, 7],
    [9, 1, 4, 1, 1, 8, 5, 6, 1],
    [1, 3, 5, 7, 2, 8, 7, 9, 1],
    [0, 9, 9, 8, 1, 6, 9, 1, 1],
    [6, 3, 2, 7, 2, 8, 0, 2, 8],
    [4, 0, 2, 3, 5, 9, 2, 5, 2]
]

spiral = spiral_order(matrix)
print("<<<" + " ".join(map(str, spiral)) + ">>>")
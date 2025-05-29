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
    [6, 5, 4, 4, 1, 2, 2, 9, 3],
    [3, 5, 5, 4, 5, 2, 8, 7, 1],
    [8, 3, 5, 5, 7, 8, 0, 4, 7],
    [7, 3, 2, 5, 3, 5, 5, 5, 8],
    [6, 4, 3, 0, 2, 2, 8, 9, 2],
    [3, 6, 4, 1, 7, 0, 6, 7, 1],
    [2, 7, 4, 3, 3, 1, 0, 6, 9],
    [7, 5, 7, 2, 2, 7, 5, 8, 2],
    [8, 8, 2, 6, 8, 9, 6, 4, 7]
]

spiral = spiral_order(matrix)
print("<<<" + " ".join(map(str, spiral)) + ">>>")
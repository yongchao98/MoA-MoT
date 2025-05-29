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
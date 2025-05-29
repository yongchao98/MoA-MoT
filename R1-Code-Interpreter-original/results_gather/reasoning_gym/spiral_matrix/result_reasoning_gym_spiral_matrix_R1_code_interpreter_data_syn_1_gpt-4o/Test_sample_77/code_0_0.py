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
    [2, 0, 3, 4, 9, 6],
    [2, 8, 2, 3, 0, 9],
    [6, 1, 8, 7, 7, 8],
    [2, 4, 7, 3, 2, 7],
    [8, 8, 8, 9, 1, 3],
    [5, 7, 8, 4, 1, 5]
]

spiral = spiral_order(matrix)
print("<<<" + " ".join(map(str, spiral)) + ">>>")
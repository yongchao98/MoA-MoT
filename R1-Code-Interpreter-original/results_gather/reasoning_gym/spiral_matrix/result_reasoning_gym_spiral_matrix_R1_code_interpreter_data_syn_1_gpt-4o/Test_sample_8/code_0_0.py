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
    [9, 5, 8, 7, 3, 0, 2, 4, 8],
    [4, 0, 3, 9, 9, 2, 4, 3, 7],
    [7, 1, 4, 3, 2, 7, 9, 0, 2],
    [2, 4, 2, 2, 7, 1, 5, 0, 6],
    [5, 9, 8, 3, 3, 9, 0, 1, 8],
    [2, 6, 2, 1, 9, 9, 3, 1, 9],
    [8, 0, 7, 8, 9, 1, 4, 9, 3],
    [5, 3, 7, 4, 8, 3, 7, 5, 1],
    [2, 4, 6, 5, 7, 7, 1, 0, 6]
]

spiral = spiral_order(matrix)
print("<<<" + " ".join(map(str, spiral)) + ">>>")
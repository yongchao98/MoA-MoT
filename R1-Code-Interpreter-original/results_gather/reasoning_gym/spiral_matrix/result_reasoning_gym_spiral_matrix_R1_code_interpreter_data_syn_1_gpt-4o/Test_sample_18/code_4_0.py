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

spiral = spiral_order(matrix)
print("<<<" + " ".join(map(str, spiral)) + ">>>")
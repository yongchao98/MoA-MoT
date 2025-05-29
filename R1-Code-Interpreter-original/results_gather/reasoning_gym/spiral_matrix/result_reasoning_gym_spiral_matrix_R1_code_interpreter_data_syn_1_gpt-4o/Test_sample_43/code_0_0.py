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
    [6, 2, 9, 3, 5, 0, 3, 8, 4],
    [6, 4, 9, 7, 4, 7, 3, 1, 4],
    [0, 4, 4, 4, 7, 6, 2, 0, 6],
    [9, 5, 7, 8, 4, 4, 7, 1, 8],
    [9, 5, 2, 9, 1, 4, 8, 7, 4],
    [8, 5, 7, 7, 6, 0, 0, 5, 2],
    [0, 3, 6, 7, 0, 5, 4, 3, 4],
    [1, 9, 8, 0, 0, 0, 2, 7, 1],
    [4, 7, 1, 7, 1, 5, 6, 8, 5]
]

spiral = spiral_order(matrix)
print("<<<" + " ".join(map(str, spiral)) + ">>>")
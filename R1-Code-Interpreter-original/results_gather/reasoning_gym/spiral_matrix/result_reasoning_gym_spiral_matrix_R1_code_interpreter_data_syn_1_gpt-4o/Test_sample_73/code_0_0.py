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
        # Take the first element of each remaining row in reverse order
        if matrix and matrix[0]:
            for row in matrix[::-1]:
                result.append(row.pop(0))
    return result

matrix = [
    [8, 7, 9, 9, 8, 5, 6],
    [2, 9, 4, 9, 8, 5, 4],
    [2, 1, 4, 5, 5, 6, 9],
    [2, 6, 0, 1, 3, 9, 5],
    [6, 4, 5, 0, 7, 2, 4],
    [6, 1, 8, 9, 7, 2, 8],
    [1, 8, 4, 2, 9, 6, 9]
]

spiral = spiral_order(matrix)
print("<<<" + " ".join(map(str, spiral)) + ">>>")
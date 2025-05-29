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
    [7, 1, 8, 7, 3, 7, 5, 5],
    [6, 1, 2, 0, 4, 8, 3, 4],
    [3, 0, 4, 5, 6, 3, 8, 1],
    [6, 7, 9, 8, 2, 7, 7, 2],
    [4, 7, 8, 7, 8, 2, 6, 3],
    [0, 8, 4, 0, 5, 9, 4, 2],
    [3, 8, 4, 9, 7, 3, 1, 5],
    [6, 6, 4, 3, 2, 2, 3, 1]
]

spiral = spiral_order(matrix)
print("<<<" + " ".join(map(str, spiral)) + ">>>")
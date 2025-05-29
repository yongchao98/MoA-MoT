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
    [6, 7, 1, 5, 2, 8],
    [6, 0, 0, 0, 3, 8],
    [5, 3, 7, 1, 3, 4],
    [6, 6, 0, 9, 8, 4],
    [2, 2, 6, 6, 5, 6],
    [2, 5, 8, 4, 3, 0]
]

spiral = spiral_order(matrix)
print("<<<" + " ".join(map(str, spiral)) + ">>>")
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
    [5, 8, 4, 3, 8, 0, 6, 4, 8],
    [3, 3, 7, 1, 1, 4, 6, 3, 9],
    [1, 9, 6, 4, 6, 2, 7, 8, 3],
    [0, 8, 2, 3, 0, 2, 1, 2, 0],
    [5, 1, 6, 1, 2, 7, 4, 6, 3],
    [7, 2, 3, 0, 0, 1, 3, 9, 3],
    [4, 3, 3, 8, 4, 2, 2, 5, 1],
    [1, 9, 0, 2, 0, 4, 0, 4, 7],
    [2, 3, 9, 0, 4, 4, 2, 0, 3]
]

spiral = spiral_order(matrix)
print("<<<" + " ".join(map(str, spiral)) + ">>>")
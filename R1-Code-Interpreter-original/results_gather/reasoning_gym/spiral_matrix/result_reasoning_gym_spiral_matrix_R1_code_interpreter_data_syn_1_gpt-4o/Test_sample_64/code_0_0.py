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
    [6, 9, 5, 3, 0, 3, 5],
    [4, 6, 8, 5, 8, 1, 4],
    [5, 7, 0, 6, 7, 3, 4],
    [1, 2, 2, 5, 5, 5, 6],
    [9, 2, 5, 7, 3, 4, 8],
    [0, 5, 4, 0, 8, 8, 9],
    [1, 8, 7, 1, 0, 1, 2]
]

spiral = spiral_order(matrix)
print("<<<" + " ".join(map(str, spiral)) + ">>>")
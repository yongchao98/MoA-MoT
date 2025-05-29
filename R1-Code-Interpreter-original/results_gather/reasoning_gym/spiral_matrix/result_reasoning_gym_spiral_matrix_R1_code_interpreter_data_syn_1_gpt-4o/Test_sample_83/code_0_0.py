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
    [1, 2, 3, 8, 1, 4],
    [3, 3, 8, 9, 8, 0],
    [3, 2, 5, 6, 5, 4],
    [4, 7, 9, 7, 5, 3],
    [9, 8, 8, 0, 2, 2],
    [0, 8, 2, 2, 4, 2]
]

spiral = spiral_order(matrix)
print("<<<" + " ".join(map(str, spiral)) + ">>>")
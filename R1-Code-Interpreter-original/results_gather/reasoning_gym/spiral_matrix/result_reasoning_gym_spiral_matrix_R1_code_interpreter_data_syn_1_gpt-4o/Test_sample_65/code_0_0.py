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
        # Take the last row, if any
        if matrix:
            result += matrix.pop()[::-1]
        # Take the first element of each remaining row, if any
        if matrix and matrix[0]:
            for row in matrix[::-1]:
                result.append(row.pop(0))
    return result

matrix = [
    [8, 0, 0, 3, 7],
    [1, 4, 5, 7, 5],
    [2, 5, 1, 6, 1],
    [3, 4, 6, 9, 8],
    [0, 0, 4, 2, 0]
]

result = spiral_order(matrix)
print("<<<" + " ".join(map(str, result)) + ">>>")
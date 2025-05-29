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
    [3, 3, 3, 1],
    [1, 4, 8, 3],
    [1, 6, 4, 5],
    [1, 2, 6, 3]
]

spiral_elements = spiral_order(matrix)
print("<<<" + " ".join(map(str, spiral_elements)) + ">>>")
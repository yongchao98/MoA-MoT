def spiral_order(matrix):
    if not matrix:
        return []

    result = []
    while matrix:
        # Add the first row
        result += matrix.pop(0)
        # Add the last column
        if matrix and matrix[0]:
            for row in matrix:
                result.append(row.pop())
        # Add the last row in reverse order
        if matrix:
            result += matrix.pop()[::-1]
        # Add the first column in reverse order
        if matrix and matrix[0]:
            for row in matrix[::-1]:
                result.append(row.pop(0))
    return result

matrix = [
    [3, 1, 3],
    [2, 4, 9],
    [1, 0, 8]
]

spiral = spiral_order(matrix)
print("<<<" + " ".join(map(str, spiral)) + ">>>")
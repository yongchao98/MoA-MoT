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
    [8, 0, 9],
    [3, 2, 9],
    [9, 8, 6]
]

spiral_list = spiral_order(matrix)
print("<<<" + " ".join(map(str, spiral_list)) + ">>>")
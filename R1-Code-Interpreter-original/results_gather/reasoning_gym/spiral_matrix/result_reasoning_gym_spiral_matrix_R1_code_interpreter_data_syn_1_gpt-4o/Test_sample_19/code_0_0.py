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
        # Add the last row
        if matrix:
            result += matrix.pop()[::-1]
        # Add the first column
        if matrix and matrix[0]:
            for row in matrix[::-1]:
                result.append(row.pop(0))
    return result

matrix = [
    [5, 6, 8, 0, 8],
    [3, 3, 6, 0, 8],
    [7, 4, 6, 2, 9],
    [0, 8, 4, 5, 4],
    [7, 1, 8, 8, 8]
]

spiral = spiral_order(matrix)
print("<<<" + " ".join(map(str, spiral)) + ">>>")
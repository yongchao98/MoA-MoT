def spiral_order(matrix):
    if not matrix:
        return []

    result = []
    while matrix:
        # Add the first row
        result += matrix.pop(0)
        # Add the last element of each remaining row
        if matrix and matrix[0]:
            for row in matrix:
                result.append(row.pop())
        # Add the last row in reverse order
        if matrix:
            result += matrix.pop()[::-1]
        # Add the first element of each remaining row (in reverse order)
        if matrix and matrix[0]:
            for row in matrix[::-1]:
                result.append(row.pop(0))
    return result

matrix = [
    [5, 0, 1, 8],
    [7, 8, 3, 6],
    [8, 5, 9, 5],
    [4, 4, 1, 1]
]

spiral_elements = spiral_order(matrix)
print("<<<" + " ".join(map(str, spiral_elements)) + ">>>")
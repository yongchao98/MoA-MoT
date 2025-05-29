def spiral_order(matrix):
    if not matrix:
        return []

    result = []
    while matrix:
        # Take the first row
        result += matrix.pop(0)
        # Take the right column
        if matrix and matrix[0]:
            for row in matrix:
                result.append(row.pop())
        # Take the bottom row
        if matrix:
            result += matrix.pop()[::-1]
        # Take the left column
        if matrix and matrix[0]:
            for row in matrix[::-1]:
                result.append(row.pop(0))
    return result

matrix = [
    [3, 2, 2, 4],
    [2, 7, 6, 7],
    [4, 8, 1, 4],
    [9, 8, 0, 2]
]

spiral_elements = spiral_order(matrix)
print("<<<" + " ".join(map(str, spiral_elements)) + ">>>")
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
    [4, 8, 7],
    [7, 4, 6],
    [4, 1, 3]
]

spiral_list = spiral_order(matrix)
print("<<<" + " ".join(map(str, spiral_list)) + ">>>")
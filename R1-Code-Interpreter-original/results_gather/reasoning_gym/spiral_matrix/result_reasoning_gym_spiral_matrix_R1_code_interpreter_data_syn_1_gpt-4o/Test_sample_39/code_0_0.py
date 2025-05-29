def spiral_order(matrix):
    if not matrix:
        return []

    result = []
    while matrix:
        # Take the first row
        result += matrix.pop(0)
        # Take the rightmost column
        if matrix and matrix[0]:
            for row in matrix:
                result.append(row.pop())
        # Take the bottom row
        if matrix:
            result += matrix.pop()[::-1]
        # Take the leftmost column
        if matrix and matrix[0]:
            for row in matrix[::-1]:
                result.append(row.pop(0))
    return result

matrix = [
    [4, 5, 3, 5, 1, 8, 9],
    [0, 9, 0, 0, 3, 2, 9],
    [3, 2, 1, 9, 9, 8, 6],
    [2, 0, 0, 4, 0, 5, 7],
    [3, 0, 8, 0, 3, 5, 8],
    [5, 4, 5, 3, 2, 7, 6],
    [9, 3, 2, 6, 8, 4, 2]
]

spiral = spiral_order(matrix)
print("<<<" + " ".join(map(str, spiral)) + ">>>")
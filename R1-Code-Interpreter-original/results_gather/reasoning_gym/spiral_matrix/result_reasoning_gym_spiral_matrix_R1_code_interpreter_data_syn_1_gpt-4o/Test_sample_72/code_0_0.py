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

# Define the matrix
matrix = [
    [3, 5],
    [6, 3]
]

# Get the spiral order
spiral = spiral_order(matrix)

# Print the result in the required format
print("<<<" + " ".join(map(str, spiral)) + ">>>")
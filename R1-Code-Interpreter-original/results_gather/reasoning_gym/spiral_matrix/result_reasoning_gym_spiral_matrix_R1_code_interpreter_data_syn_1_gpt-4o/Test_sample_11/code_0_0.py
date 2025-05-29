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

# Define the matrix
matrix = [
    [5, 7],
    [2, 4]
]

# Get the spiral order
spiral = spiral_order(matrix)

# Print the result in the required format
print("<<<" + " ".join(map(str, spiral)) + ">>>")
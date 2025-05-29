def spiral_order(matrix):
    if not matrix or not matrix[0]:
        return []

    top, bottom = 0, len(matrix) - 1
    left, right = 0, len(matrix[0]) - 1
    result = []

    while top <= bottom and left <= right:
        # Traverse from left to right along the top row
        for i in range(left, right + 1):
            result.append(matrix[top][i])
        top += 1

        # Traverse from top to bottom along the right column
        for i in range(top, bottom + 1):
            result.append(matrix[i][right])
        right -= 1

        if top <= bottom:
            # Traverse from right to left along the bottom row
            for i in range(right, left - 1, -1):
                result.append(matrix[bottom][i])
            bottom -= 1

        if left <= right:
            # Traverse from bottom to top along the left column
            for i in range(bottom, top - 1, -1):
                result.append(matrix[i][left])
            left += 1

    return result

# Define the matrix
matrix = [
    [7, 8, 3, 7],
    [8, 5, 1, 6],
    [7, 1, 2, 6],
    [4, 0, 5, 1]
]

# Get the spiral order
spiral_list = spiral_order(matrix)

# Print the result in the required format
print("<<<" + " ".join(map(str, spiral_list)) + ">>>")
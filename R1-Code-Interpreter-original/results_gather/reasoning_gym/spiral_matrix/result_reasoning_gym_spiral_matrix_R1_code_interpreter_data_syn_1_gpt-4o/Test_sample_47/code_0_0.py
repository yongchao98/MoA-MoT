def spiral_order(matrix):
    if not matrix:
        return []
    
    top, bottom = 0, len(matrix) - 1
    left, right = 0, len(matrix[0]) - 1
    result = []
    
    while top <= bottom and left <= right:
        # Traverse from left to right across the top row
        for i in range(left, right + 1):
            result.append(matrix[top][i])
        top += 1
        
        # Traverse from top to bottom down the right column
        for i in range(top, bottom + 1):
            result.append(matrix[i][right])
        right -= 1
        
        if top <= bottom:
            # Traverse from right to left across the bottom row
            for i in range(right, left - 1, -1):
                result.append(matrix[bottom][i])
            bottom -= 1
        
        if left <= right:
            # Traverse from bottom to top up the left column
            for i in range(bottom, top - 1, -1):
                result.append(matrix[i][left])
            left += 1
    
    return result

# Given matrix
matrix = [
    [3, 4, 1, 1],
    [7, 2, 9, 0],
    [7, 5, 5, 3],
    [3, 2, 2, 7]
]

# Get the spiral order
spiral_list = spiral_order(matrix)

# Format the output
print("<<<" + " ".join(map(str, spiral_list)) + ">>>")
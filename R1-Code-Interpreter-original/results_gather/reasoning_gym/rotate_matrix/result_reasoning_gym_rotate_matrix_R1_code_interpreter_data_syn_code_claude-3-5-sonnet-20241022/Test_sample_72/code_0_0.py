def rotate_matrix(matrix, degrees):
    # Convert degrees to number of 90-degree rotations
    rotations = (degrees // 90) % 4
    
    # If rotations is 0, return the original matrix
    if rotations == 0:
        # Format the output as required
        rows = len(matrix)
        cols = len(matrix[0])
        result = ""
        for i in range(rows):
            for j in range(cols):
                result += str(matrix[i][j])
                if j < cols - 1:
                    result += " "
            if i < rows - 1:
                result += "\n"
        print(result)

# Test with the given matrix
matrix = [[3, 3], [6, 5]]
rotate_matrix(matrix, 0)
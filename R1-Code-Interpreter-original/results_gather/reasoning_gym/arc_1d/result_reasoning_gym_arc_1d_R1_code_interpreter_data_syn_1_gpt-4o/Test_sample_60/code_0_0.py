def transform_grid(input_grid):
    # Convert the input list to a 5x5 matrix
    matrix = [input_grid[i:i+5] for i in range(0, 25, 5)]
    
    # Initialize the output matrix with zeros
    output_matrix = [[0]*5 for _ in range(5)]
    
    # Copy the first row directly
    output_matrix[0] = matrix[0]
    
    # Fill the rest of the matrix based on the rule
    for i in range(1, 5):
        for j in range(5):
            if matrix[i][j] != 0:
                # Fill the entire row and column with the non-zero number
                for k in range(5):
                    output_matrix[i][k] = matrix[i][j]
                    output_matrix[k][j] = matrix[i][j]
    
    # Flatten the output matrix to a list
    output_grid = [num for row in output_matrix for num in row]
    return output_grid

# Test input
test_input = [3, 3, 3, 0, 0, 0, 3, 0, 0, 0, 3, 0, 0, 0, 3, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0]
output = transform_grid(test_input)
print(output)
